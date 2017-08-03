#include <rcppExport.h>
#include <combinatorics.h>
#include <dataStructure.h>
#include <types.h>
#include <zdensity.h>
#include <bfgs.h>
#include <optimize.h>
#include <fpUcHandling.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <iostream>
#include <climits>
#include <cmath> // for the long double overloading of exp and log
#include <sstream>
#include <string>
#include <stdexcept>

// using pretty much:
using std::map;
using std::set;
using std::vector;
using std::accumulate;
using std::find;
using std::set_difference;
using std::count;
using std::max_element;

using namespace Rcpp;

// ***************************************************************************************************//

// compute varying part of log marginal likelihood for specific GLM / Cox model
// plus byproducts.
double
getGlmVarLogMargLik(const ModelPar &mod,
                    const DataValues& data,
                    const FpInfo& fpInfo,
                    const UcInfo& ucInfo,
                    const FixInfo& fixInfo,
                    const Book& bookkeep,
                    const GlmModelConfig& config,
                    const GaussHermite& gaussHermite,
                    Cache& cache,
                    double& zMode,
                    double& zVar,
                    double& laplaceApprox,
                    double& residualDeviance)
{
    // echo detailed progress in debug mode
    if(bookkeep.debug)
    {
        Rprintf("\ngetGlmVarLogMargLik: Starting log marginal likelihood approximation for:\n%s", mod.print(fpInfo).c_str());
    }

    // check if any interrupt signals have been entered
    R_CheckUserInterrupt();

    // the return value will be placed in here:
    double ret = 0.0;

    // protect everything for problems in IWLS
    try
    {
        // get negative log unnormalized z density: a function object.
        NegLogUnnormZDens negLogUnnormZDens(mod,
                                            data,
                                            fpInfo,
                                            ucInfo,
                                            fixInfo,
                                            config,
                                            bookkeep);
        residualDeviance = negLogUnnormZDens.getResidualDeviance();

        // try to ask for analytic solutions in the TBF case
        // (they are available for the incomplete inverse gamma hyperprior and of course
        // for the fixed g case)
        if(bookkeep.tbf)
        {
            if(bookkeep.debug)
            {
                Rprintf("\ngetGlmVarLogMargLik: TBF used, ask for analytic solution ...");
            }

            if(bookkeep.useFixedg)
            {
                zMode = log(config.fixedg);
                ret = - negLogUnnormZDens(zMode);

                if(bookkeep.debug)
                {
                    Rprintf("\ngetGlmVarLogMargLik: using fixed value g = %f", config.fixedg);
                }
            }
            else if(bookkeep.empiricalBayes)
            {
                ret = negLogUnnormZDens.getTBFMaxLogCondMargLik(zMode);

                if(bookkeep.debug)
                {
                    Rprintf("\ngetGlmVarLogMargLik: local empirical Bayes estimate of g is %f", exp(zMode));
                }
            }
            else
            {
                ret = negLogUnnormZDens.getTBFLogMargLik();
            }

            if(! R_IsNA(ret))
            {
                if(bookkeep.debug)
                {
                    Rprintf("\ngetGlmVarLogMargLik: analytic solution was found with result %f", ret);
                }
                return ret;
            }
            else
            {
                if(bookkeep.debug)
                {
                    Rprintf("\ngetGlmVarLogMargLik: no analytic solution possible for this hyperprior on g");
                }
            }

            // if no analytic solutions are available, we continue with the code below
        }

        // cache this function, because we do not want to evaluate it more often
        // than necessary.
        CachedFunction<NegLogUnnormZDens> cachedNegLogUnnormZDens(negLogUnnormZDens);

        if(bookkeep.useFixedg)
        {
            zMode = log(config.fixedg);

            // echo detailed progress in debug mode
            if(bookkeep.debug)
            {
                Rprintf("\ngetGlmVarLogMargLik: returning conditional marginal likelihood at fixed z = %f", zMode);
            }

            // return the log conditional marginal density log f(y | zfixed)
            ret = - cachedNegLogUnnormZDens(zMode);
        }
        else if(bookkeep.empiricalBayes)
        {
            // construct an appropriate object for using the optimize routine
            Brent<CachedFunction<NegLogUnnormZDens> > brent(cachedNegLogUnnormZDens,
                                                            -20,//-100.0, // log(DBL_MIN) + 40.0,
                                                            20,//+200.0, // log(DBL_MAX) - 40.0,
                                                            sqrt(EPS));
          
          // for(int zz = -55; zz < 40; zz=zz+2){
          //   Rcpp::Rcout << zz<< " : "<< cachedNegLogUnnormZDens(zz) << std::endl;
          // }
            // and get the mode from that.
            zMode = brent.minimize();

          // Rcpp::Rcout <<  "Zmode is : "<< zMode << std::endl;
          
            // zVar and laplaceApprox are not touched, as they are not needed!

            // echo detailed progress in debug mode
            if(bookkeep.debug)
            {
                Rprintf("\ngetGlmVarLogMargLik: finished optimization of conditional marginal likelihood at mode %f", zMode);
            }

            // return the log conditional marginal density log f(y | z_mode)
            ret = - cachedNegLogUnnormZDens(zMode);
        }
        else // start full Bayes
        {
            // get function invHess to compute an accurate variance estimate
            AccurateNumericInvHessian<CachedFunction<NegLogUnnormZDens> > invHess(cachedNegLogUnnormZDens);

            // decide if bfgs, or optimize should be used
            if(bookkeep.useBfgs)
            {
                // and run the minimization algorithm on it.

                // constrain z to lie in the interval [log(DBL_MIN), log(DBL_MAX)], so
                // that g = exp(z) is always in [DBL_MIN, DBL_MAX].
                // and we put a little safety margin on it.

                // 18/05: test a bit more realistic interval (+-200)
                // 08/07: even further shorten the interval to -100, 200 because the posterior
                //        will be very flat below z=-100
                Bfgs<CachedFunction<NegLogUnnormZDens> > bfgs(cachedNegLogUnnormZDens,
                                                              bookkeep.debug,
                                                              -100.0, // log(DBL_MIN) + 40.0,
                                                              +200.0); //log(DBL_MAX) - 40.0);

                // now minimize the negative log density,
                // and put the resulting mode into zMode.
                int convergence = bfgs.minimize(0.0,
                                                zMode,
                                                zVar);

                // if we lost precision, do a second minimization starting from the previous mode
                // and estimate the variance afterwards separately (otherwise the variance estimate
                // would not be reliable!)
                if(convergence == -1)
                {
                    bfgs.minimize(zMode,
                                  zMode,
                                  zVar);
                    zVar = invHess(zMode);
                }
            }
            else // use optimize
            {
                // construct an appropriate object for using the optimize routine
                Brent<CachedFunction<NegLogUnnormZDens> > brent(cachedNegLogUnnormZDens,
                                                                -100.0, // log(DBL_MIN) + 40.0,
                                                                +200.0, // log(DBL_MAX) - 40.0,
                                                                sqrt(EPS));
                // and get the mode from that.
                zMode = brent.minimize();

                // here we have to compute the inverse hessian afterwards:
                // use the same epsilon here as for the minimization routine.
                zVar = invHess(zMode);
            }

            // be careful that the result for zVar is not totally wrong.
            // if the variance estimate is very large, warn the user.
            if(zVar > bookkeep.largeVariance)
            {
                Rf_warning("\nLarge variance estimate (%f > %f) for z marginal encountered by BFGS",
                           zVar, bookkeep.largeVariance);
            }

            // since we minimize here the negative log density, the resulting inverse
            // hessian is the estimated variance of the approximating Gaussian.

            // be sure that we get a positive variance estimate
            if(zVar <= 0.0)
            {
                // warn
                Rf_warning("\nNon-positive variance estimate %f encountered!\nProbably the optimization did not converge.\nResetting to default variance",
                           zVar);

                // set large default to explore some space
                zVar = 15.0;
            }

            // echo detailed progress in debug mode
            if(bookkeep.debug)
            {
                Rprintf("\ngetGlmVarLogMargLik: finished minimization at mode %f with variance %f", zMode, zVar);
            }

            // compute the Laplace approximation to the log marginal likelihood from the mode and variance
            laplaceApprox = M_LN_SQRT_2PI + 0.5 * log(zVar) - cachedNegLogUnnormZDens(zMode);
            // so this does not require evaluations inside the Gauss-Hermite quadrature.

            // then compute the Gauss-Hermite quadrature, using the supplied standard nodes and
            // weights from R
            MyDoubleVector nodes;
            MyDoubleVector logWeights;

            // get the nodes and log weights for this mode and variance:
            gaussHermite.getNodesAndLogWeights(zMode, zVar, nodes, logWeights);

            // the log contributions which will be stored here:
            SafeSum logContributions;

            // compute them now
            MyDoubleVector::const_iterator n = nodes.begin();
            for(MyDoubleVector::const_iterator
                    w = logWeights.begin();
                    w != logWeights.end();
                    ++w, ++n)
            {
                logContributions.add((*w) - cachedNegLogUnnormZDens(*n));
            }

            // the result is the log of the sum of the exp'ed values
            ret = logContributions.logSumExp();
        } // end full Bayes

        // echo detailed progress in debug mode
        if(bookkeep.debug)
        {
            Rprintf("\ngetGlmVarLogMargLik: finished log marginal likelihood approximation.");
        }

        // also give back the cache
        cache = cachedNegLogUnnormZDens.getCache();

        // check finiteness
        if (! R_finite(ret))
        {
            std::ostringstream stream;
            stream << "getGlmVarLogMargLik got non-finite log marginal likelihood approximation " << ret;
            throw std::domain_error(stream.str().c_str());
        }
    }
    catch (std::domain_error& error)
    {
        if(bookkeep.debug)
        {
            Rprintf("\ngetGlmVarLogMargLik: Model can not be included, because\n%s", error.what());
        }

        return R_NaN;
    }

    // only then return the estimate
    return ret;
}

// ***************************************************************************************************//

// compute varying part of logarithm of model prior
double
getVarLogPrior(
               const ModelPar &mod,
               const FpInfo &fpInfo,
               const UcInfo &ucInfo,
               const FixInfo& fixInfo,
               const Book &bookkeep)
{
    if (bookkeep.modelPrior == "sparse")
    {
        SafeSum thisVarLogPrior;
        for (unsigned int i = 0; i != fpInfo.nFps; i++)
        { // for each fp covariate
            unsigned int degree = mod.fpPars.at(i).size();
            double thisVal = -Rf_lchoose(fpInfo.fpcards[i] - 1 + degree, degree)
                    - log1p(fpInfo.fpmaxs[i]);
            thisVarLogPrior.add(thisVal);
        }
        return thisVarLogPrior.sum() - ((ucInfo.nUcGroups + fixInfo.nFixGroups) * M_LN2); //TODO Is this correct?
    }
    else if (bookkeep.modelPrior == "dependent")
    {
        // determine number of all covariates (covariate groups):
        int nCovs = ucInfo.nUcGroups + fixInfo.nFixGroups + fpInfo.nFps; //TODO Is this correct?

        // determine number of included FPs and which are nonlinear:
        int nInclContinuous = 0;
        PosIntVector nonlinearFps;

        for(PosInt i = 0; i != fpInfo.nFps; i++)
        {
            Powers powersi = mod.fpPars.at(i);
            if (! powersi.empty())
            {
                ++nInclContinuous;

                if(mod.fpPars.at(i) != fpInfo.linearPowers)
                {
                    nonlinearFps.push_back(i);
                }
            }
        }

        // determine number of included discrete covariates:
        int nInclDiscrete = mod.ucPars.size() + fixInfo.nFixGroups;

        // so altogether there are
        int nIncluded = nInclContinuous + nInclDiscrete;
        // included covariates

        // and the number of possible nonlinear transformations
        // for each variable is also included in the computations:
        double sumLogNonlinearPossibilities = 0.0;
        for(PosIntVector::const_iterator
                i = nonlinearFps.begin();
                i != nonlinearFps.end();
                ++i)
        {
            sumLogNonlinearPossibilities += log(fpInfo.numberPossibleFps.at(*i) - 2.0);
            //                              Note: degree 0 and linear degree 1 FP are subtracted
        }

        double result = - log1p(nCovs) - Rf_lchoose(nCovs, nIncluded) -
                log1p(nInclContinuous) - Rf_lchoose(nInclContinuous, nonlinearFps.size()) -
                sumLogNonlinearPossibilities;

        return result;
    }
    else
    {
        return - (ucInfo.nUcGroups * M_LN2);
    }
}

// ***************************************************************************************************//

// 21/11/2012: modify for tbf methodology

// compute (varying part of) marginal likelihood and prior of model and insert it into model set
void
computeGlm(const ModelPar &mod,
           set<Model> &space,
           const DataValues& data,
           const FpInfo& fpInfo,
           const UcInfo& ucInfo,
           const FixInfo& fixInfo,
           Book& bookkeep,
           const GlmModelConfig& config,
           const GaussHermite& gaussHermite)
{
    // log prior
    const double thisLogPrior = getVarLogPrior(mod,
                                               fpInfo,
                                               ucInfo,
                                               fixInfo,
                                               bookkeep);

    // initialize variables which will definitely be overwritten below
    double thisVarLogMargLik = 0.0;

    // this will not be available when this is the null model! then it is *not* overwritten below.
    double zMode = R_NaReal;
    double zVar = R_NaReal;
    double laplaceApprox = R_NaReal;
    double residualDeviance = R_NaReal;
    Cache cache;

    // compute log marginal likelihood, and also as byproducts unnormalized z density information.

    // be careful: if this is the null model, then just input the data computed in R.
    if(mod.size(ucInfo, fixInfo) == 0)
    {
        thisVarLogMargLik = config.nullModelLogMargLik;
    }
    else // not the null model, so at least one other coefficient than the intercept present in the model
    {
        thisVarLogMargLik = getGlmVarLogMargLik(mod, data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite,
                                                cache, zMode, zVar, laplaceApprox, residualDeviance);
    }

    // if we get back NaN
    if (R_IsNaN(thisVarLogMargLik) == TRUE)
    {
        // increment counter of bad models
        bookkeep.nanCounter++;
        // we do not save this here, because for the sampling mode we will have a different function!
    }
    else
    {
        // put all into the modelInfo
        GlmModelInfo info(thisVarLogMargLik, thisLogPrior, cache, zMode, zVar, laplaceApprox, residualDeviance);

        // altogether we have the model:
        Model thisModel(mod, info);

        // and insert it into the model space
        if (space.size() >= bookkeep.nModels)
        {
            set<Model>::iterator it = space.begin();
            if (* it < thisModel)
            { // compare this model to the least probable model in the set
                space.erase(it);
                space.insert(thisModel); // exchange if it is better than this worst model in the set
            }
        }
        else
        {
            space.insert(thisModel);
        }

        // compute log posterior probability (up to an additive constant) and
        // append it to the safe sum object
        bookkeep.modelLogPosteriors.add(thisVarLogMargLik + thisLogPrior);

        // update inclusion probabilities for covariate (groups)
        mod.pushInclusionProbs(fpInfo, ucInfo, bookkeep); //TODO

        // increment distinct models counter
        bookkeep.modelCounter++;
    }

    // display computation progress at each percent:
    if (((bookkeep.modelCounter + bookkeep.nanCounter) %
            std::max(data.totalNumber / 100,
                static_cast<PosLargeInt> (1))
                == 0)
            && bookkeep.verbose) // but only print if verbose option was chosen
    {
        Rprintf("-");
    }
}


// ***************************************************************************************************//


//TODO FINISH THIS FUNCTION

// compute only the models in the list R-list "R_modelConfigs"
List
glmModelsInList(const DataValues& data,
                const FpInfo& fpInfo,
                const UcInfo& ucInfo,
                const FixInfo& fixInfo,
                Book& bookkeep,
                const GlmModelConfig& config,
                const GaussHermite& gaussHermite,
                const List& rcpp_modelConfigs)
{
    // ------------
    // bookkeeping:

    // for computation of inclusion probs:
    // array of IndexSafeSum objects.
//     const int cgwp_length = fpInfo.nFps + ucInfo.nUcGroups;
//     //IndexSafeSum cgwp[fpInfo.nFps + ucInfo.nUcGroups];
//     IndexSafeSum cgwp[cgwp_length];
//     bookkeep.covGroupWisePosteriors = cgwp;

    std::vector<IndexSafeSum> cgwp(fpInfo.nFps + ucInfo.nUcGroups);
    bookkeep.covGroupWisePosteriors = cgwp;
    
    
    
    // start the set of ordered models:
    set<Model> orderedModels;

    // check that the set is large enough for all models in the list
    if (orderedModels.max_size() < static_cast<set<Model>::size_type >(rcpp_modelConfigs.size()))
        Rf_error("\nlist of model space is too large - cannot compute every model\n");

    // ------------
    // process all models in the list:

    for(R_len_t i = 0; i < rcpp_modelConfigs.size(); ++i)
    {
    as<List>(rcpp_modelConfigs[i]);
    
        // this is the current model config:
        ModelPar modelConfig(as<List>(rcpp_modelConfigs[i]),
                             fpInfo);

        // compute this one
        computeGlm(modelConfig, orderedModels,
                   data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);

    }

    // ------------
    // we have finished.

    // now echo statistics:
    if (bookkeep.verbose)
    {
        Rprintf("\nActual number of possible models:  %d ",
                bookkeep.modelCounter);
        Rprintf("\nNumber of non-identifiable models: %d",
                bookkeep.nanCounter);
        Rprintf("\nNumber of saved possible models:   %d\n",
                orderedModels.size());
    }

    // ------------
    // allocate the return list
    List ret(orderedModels.size());

    // and fill it:

    // normalize posterior probabilities and correct log prior of the models to return
    // (we do not know here the normalizing constant for the marginal likelihoods!)
    const long double logNormConst = bookkeep.modelLogPosteriors.logSumExp();

    // first the single models
    PosInt i = 0;
    for (set<Model>::const_reverse_iterator
            j = orderedModels.rbegin();
            j != orderedModels.rend();
            j++)
    {
        ret[i++] = j->convert2list(fpInfo,
                                   logNormConst,
                                   bookkeep);
    }

    // then some attributes:

    NumericVector inc(fpInfo.nFps + ucInfo.nUcGroups);
    for (R_len_t i = 0; i != inc.size(); ++i)
    {
        inc[i] = bookkeep.covGroupWisePosteriors[i].sumNormalizedExp(bookkeep.modelLogPosteriors, logNormConst);
    }
    ret.attr("inclusionProbs") = inc;
    ret.attr("numVisited") = static_cast<double>(bookkeep.modelCounter);
    ret.attr("logNormConst") = static_cast<double>(logNormConst);

    // ------------
    // finally return the list.
    return ret;
}


// ***************************************************************************************************//

// recursion via:
void
glmPermPars(PosInt pos, // current position in parameter vector, starting from 0 - copied.
            ModelPar mod, // is copied every time! everything else is call by reference.
            set<Model>& space, // the model space
            const DataValues& data,
            const FpInfo& fpInfo,
            const UcInfo& ucInfo,
            const FixInfo& fixInfo,
            Book& bookkeep,
            const GlmModelConfig& config,
            const GaussHermite& gaussHermite)

{
    // if some fps are still left
    if (pos != fpInfo.nFps)
    {
        // cardinality of the power set at this position
        const PosInt card = fpInfo.fpcards.at(pos);

        // degree 0:
        glmPermPars(pos + 1, mod, space,
                    data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);

        // different degrees for fp at pos:
        // degrees 1, ..., fpmax
        for (PosInt deg = 1; deg <= fpInfo.fpmaxs.at(pos); deg++)
        {
            // increment sums of fp degrees
            mod.fpSize++;

            // partition of deg into card parts
            IntVector part (card);

            // internal variables for comp_next
            bool more1 = false;
            int h(0), t(0);

            do
            {
                // next partition of deg into card parts
                comp_next(deg, card, part, &more1, h, t);

                // convert into multiset
                mod.fpPars.at(pos) = freqvec2Powers(part, card);

                // and go on
                glmPermPars(pos + 1, mod, space,
                            data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);
            }
            while (more1);
        }
    }
    else // no fps left (all FPs have received their powers)
    {
        // no uc group
        computeGlm(mod, space,
                   data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite); //TODO IS THIS THE NULL  MODEL? IF SO ADD NULL+FIXED NEXT

        // different positive number (deg) of uc groups
        for (PosInt deg = 1; deg <= ucInfo.nUcGroups; deg++)
        {
            // partition of deg into card parts
          IntVector subset (deg);

            // internal variables for ksub_next
            bool more2 = false;
            int m(0), m2(0);

            do
            {
                // next subset (positive integers)
                ksub_next(ucInfo.nUcGroups, deg, subset, &more2, m, m2);

                // convert into set
                mod.ucPars = IntSet(subset.data(), subset.data() + deg);

                // and compute this model
                computeGlm(mod, space,
                           data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);
            }
            while (more2);
        }
    }
}

// ***************************************************************************************************//

List
glmSampling(const DataValues& data,
            const FpInfo& fpInfo,
            const UcInfo& ucInfo,
            const FixInfo& fixInfo,
            Book& bookkeep,
            const GlmModelConfig& config,
            const GaussHermite& gaussHermite)
{
    // models which can be found during chain run can be cached in here:
    ModelCache modelCache(bookkeep.nCache);

    // upper limit for num of columns: min(n, maximum fixed + fp + uc columns).
    PosInt maxDim = std::min(static_cast<PosInt>(data.nObs), 1 + fpInfo.maxFpDim + ucInfo.maxUcDim);

    // the FP range
    const PosIntSet fpRange = constructSequence(fpInfo.nFps);

    // start model is the null model!
    ModelMcmc old(fpInfo,
                  ucInfo,
                  maxDim,
                  config.nullModelLogMargLik);

    // insert this model into the cache
    double logPrior = getVarLogPrior(old.modPar,
                                     fpInfo,
                                     ucInfo,
                                     fixInfo,
                                     bookkeep);
    old.logPrior = logPrior;

    // put all into the modelInfo
    GlmModelInfo startInfo(old.logMargLik, logPrior, Cache(), R_NaReal, R_NaReal, R_NaReal, 0.0);

    modelCache.insert(old.modPar, startInfo);

    // start with this model config
    ModelMcmc now(old);

    if(fixInfo.nFixGroups > 0){
      // move to the null model + fixed covariates **********************************************//
      
      // add the fixed covariates to the model configuration
      IntSet s;
      for (unsigned int i = 0; i < fixInfo.nFixGroups; ++i) 
        s.insert(s.end(), i+1);
      now.modPar.fixPars = s;
      
      //get log prior for null+fixed model
     double logPrior2 = getVarLogPrior(now.modPar,
                                fpInfo,
                                ucInfo,
                                fixInfo,
                                bookkeep);
      
      // and marginal log like
      double zMode = 0.0;
      double zVar = 0.0;
      double laplaceApprox = 0.0;
      double residualDeviance = R_NaReal;
      Cache cache;
      
      now.logMargLik = getGlmVarLogMargLik(now.modPar,
                                           data,
                                           fpInfo,
                                           ucInfo,
                                           fixInfo,
                                           bookkeep,
                                           config,
                                           gaussHermite,
                                           cache,
                                           zMode,
                                           zVar,
                                           laplaceApprox,
                                           residualDeviance);
      
      // put all into the modelInfo
      GlmModelInfo start2Info(now.logMargLik, logPrior2, Cache(), R_NaReal, R_NaReal, R_NaReal, 0.0);
      
      modelCache.insert(now.modPar, start2Info);
      
      // start with this model config
      ModelMcmc now2(now);
      
      now = now2;
    }
    
    // Start MCMC sampler***********************************************************//

    GetRNGstate(); // use R's random number generator

    for(PosLargeInt t = 0; t != bookkeep.chainlength; /* ++t explicitly at the end */)
    {
            double logPropRatio; // log proposal ratio

            // randomly select move type
            double u1 = unif_rand();

            if (u1 < old.birthprob)
            {                                                                                        // BIRTH
                    PosInt newCovInd = *discreteUniform<PosIntSet>(old.freeCovs);

                    if (newCovInd <= fpInfo.nFps)
                    {                                                                                // some fp index
                            Int powerIndex = discreteUniform<Int>(0, fpInfo.fpcards[newCovInd-1]);
                            now.modPar.fpPars.at(newCovInd-1).insert(powerIndex);
                            now.modPar.fpSize++; // correct invariants
                            now.dim++;
                            PosInt newPowersEqualPowerIndex = count(now.modPar.fpPars.at(newCovInd-1).begin(), now.modPar.fpPars.at(newCovInd-1).end(), powerIndex);
                            PosInt m = old.modPar.fpPars.at(newCovInd-1).size();

                            logPropRatio = log(double(newPowersEqualPowerIndex)) + log(double(fpInfo.fpcards[newCovInd-1])) - log1p(m);
                    }
                    else
                    {                                                                                // uc index
                            Int index = *discreteUniform(old.freeUcs);
                            now.modPar.ucPars.insert(index);
                            now.dim += ucInfo.ucSizes.at(index - 1);
                            now.freeUcs = now.modPar.getFreeUcs(ucInfo.ucSizes, now.dim, maxDim);

                            logPropRatio = log(double(old.freeUcs.size())) - log(double(now.modPar.ucPars.size()));
                    }

                    now.presentCovs.insert(newCovInd);
                    now.freeCovs = now.modPar.getFreeCovs(fpInfo, now.freeUcs, now.dim, maxDim);

                    if (now.dim == maxDim)
                    {
                            now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3 : 0.5;
                    } else
                    {
                            now.birthprob = now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3;
                    }

                    logPropRatio += log(now.deathprob) - log(old.birthprob) + log(double(old.freeCovs.size())) - log(double(now.presentCovs.size()));

            }
            else if
            (u1 < old.birthprob + old.deathprob)
            {                                                                                      // DEATH
                    PosInt oldCovInd = *discreteUniform(old.presentCovs);

                    if (oldCovInd <= fpInfo.nFps)
                    {                                                                            // some fp index
                            Powers::iterator powerIterator = discreteUniform(now.modPar.fpPars.at(oldCovInd-1));
                            PosInt oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(oldCovInd-1).begin(), old.modPar.fpPars.at(oldCovInd-1).end(), *powerIterator);
                            now.modPar.fpPars.at(oldCovInd-1).erase(powerIterator);
                            now.modPar.fpSize--; // correct invariants
                            now.dim--;

                            logPropRatio =  - log(double(oldPowersEqualPowerIndex)) - log(double(fpInfo.fpcards[oldCovInd-1])) + log(double(old.modPar.fpPars.at(oldCovInd-1).size()));

                    } else {                                                                                                        // uc index
                            IntSet::iterator IndIterator = discreteUniform(now.modPar.ucPars);
//                            now.modPar.ucSize--;
                            now.dim -= ucInfo.ucSizes.at(*IndIterator - 1);
                            now.modPar.ucPars.erase(IndIterator);
                            now.freeUcs = now.modPar.getFreeUcs(ucInfo.ucSizes, now.dim, maxDim);
                            logPropRatio = log(double(old.modPar.ucPars.size())) - log(double(now.freeUcs.size()));
                    }
                    now.presentCovs = now.modPar.getPresentCovs();
                    now.freeCovs = now.modPar.getFreeCovs(fpInfo, now.freeUcs, now.dim, maxDim);
                    if (now.dim == 1)
                    {
                            now.birthprob = 1; now.deathprob = now.moveprob = 0;
                    } else
                    {
                            now.birthprob = now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3.0;
                    }
                    logPropRatio += log(now.birthprob) - log(old.deathprob) + log(double(old.presentCovs.size())) - log(double(now.freeCovs.size()));

            }
            else if (u1 < old.birthprob + old.deathprob + old.moveprob)
            {                                                                                   // MOVE
                    PosInt CovInd = *discreteUniform<PosIntSet>(old.presentCovs);

                    if (CovInd <= fpInfo.nFps)
                    {                                                                     // some fp index
                            Powers::iterator powerIterator = discreteUniform(now.modPar.fpPars.at(CovInd-1));
                            PosInt oldPowersEqualPowerIndex = count(old.modPar.fpPars.at(CovInd-1).begin(), old.modPar.fpPars.at(CovInd-1).end(), *powerIterator);
                            now.modPar.fpPars.at(CovInd-1).erase(powerIterator);
                            Int powerIndex = discreteUniform<Int>(0, fpInfo.fpcards[CovInd-1]);
                            now.modPar.fpPars.at(CovInd-1).insert(powerIndex);
                            PosInt newPowersEqualPowerIndex = count(now.modPar.fpPars.at(CovInd-1).begin(), now.modPar.fpPars.at(CovInd-1).end(), powerIndex);
                            // free, present Covs and move type probs are unchanged
                            logPropRatio = log(double(newPowersEqualPowerIndex)) - log(double(oldPowersEqualPowerIndex));
                    }
                    else
                    {                                                                                                        // uc index
                            IntSet::iterator IndIterator = discreteUniform(now.modPar.ucPars);
                            now.dim -= ucInfo.ucSizes.at(*IndIterator - 1);
                            now.modPar.ucPars.erase(IndIterator);
                            now.freeUcs = now.modPar.getFreeUcs(ucInfo.ucSizes, now.dim, maxDim);
                            Int index = *discreteUniform<IntSet>(now.freeUcs);
                            now.modPar.ucPars.insert(index);
                            now.dim += ucInfo.ucSizes.at(index - 1);
                            now.freeUcs = now.modPar.getFreeUcs(ucInfo.ucSizes, now.dim, maxDim);
                            // here something may change, therefore:
                            now.freeCovs = now.modPar.getFreeCovs(fpInfo, now.freeUcs, now.dim, maxDim);
                            if (now.dim == maxDim)
                            {
                                    now.birthprob = 0; now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 1.0 / 3.0 : 0.5;
                            }
                            else
                            {
                                    now.birthprob = now.deathprob = now.moveprob = (now.modPar.fpSize > 0) ? 0.25 : 1.0 / 3.0;
                            }
                            logPropRatio = 0.0;
                    }
            } else {                                                                                                        // SWITCH (of FP vectors)
                    // select only the FP present covs
                    PosIntSet presentFps = removeElement(old.presentCovs, fpInfo.nFps + 1);

                    // so we have the first power vector:
                    PosInt firstFpInd = *discreteUniform<PosIntSet>(presentFps);
                    Powers first = now.modPar.fpPars.at(firstFpInd - 1);

                    // the second power vector from all other FPs
                    PosIntSet otherFps = removeElement(fpRange, firstFpInd);
                    PosInt secondFpInd = *discreteUniform<PosIntSet>(otherFps);
                    Powers second = now.modPar.fpPars.at(secondFpInd - 1);

                    // save the first
                    Powers saveFirst = first;

                    // copy second to first
                    now.modPar.fpPars.at(firstFpInd - 1) = second;

                    // and save to second
                    now.modPar.fpPars.at(secondFpInd - 1) = saveFirst;

                    // so now we have switched the power vectors.

                    // move type probs are not changed, because the number of present FPs is unchanged,
                    // as well as the dimension of the model.

                    // but carefully update the information which covariates are free and which are present
                    now.freeCovs = now.modPar.getFreeCovs(fpInfo, now.freeUcs, now.dim, maxDim);
                    now.presentCovs = now.modPar.getPresentCovs();

                    // and the proposal ratio is 1, thus the log proposal ratio is 0:
                    logPropRatio = 0;
            }

            // search for log marg lik of proposed model
            GlmModelInfo nowInfo = modelCache.getModelInfo(now.modPar);

            if (R_IsNA(nowInfo.logMargLik))
            { // "now" is a new model

                double zMode = 0.0;
                double zVar = 0.0;
                double laplaceApprox = 0.0;
                double residualDeviance = R_NaReal;
                Cache cache;

                // so we must compute the log marg lik now.
                now.logMargLik = getGlmVarLogMargLik(now.modPar,
                                                 data,
                                                 fpInfo,
                                                 ucInfo,
                                                 fixInfo,
                                                 bookkeep,
                                                 config,
                                                 gaussHermite,
                                                 cache,
                                                 zMode,
                                                 zVar,
                                                 laplaceApprox,
                                                 residualDeviance);

                // check if the new model is OK
                if (R_IsNaN(now.logMargLik))
                {
                    // we do not save this model in the model cache
                    bookkeep.nanCounter++;
                }
                else
                { // OK: then compute the rest, and insert into model cache

                    now.logPrior = getVarLogPrior(now.modPar,
                                              fpInfo,
                                              ucInfo,
                                              fixInfo,
                                              bookkeep);

                    // insert the model parameter/info into the model cache

                    // problem: this could erase the old model from the model cache,
                    // and invalidate the iterator old.mapPos!
                    // ==> so we cannot work with the iterators here.
                    modelCache.insert(now.modPar,
                                  GlmModelInfo(now.logMargLik,
                                               now.logPrior,
                                               cache,
                                               zMode,
                                               zVar,
                                               laplaceApprox,
                                               residualDeviance));
                }
            }
            else // "now" is an old model
            {
                // extract log marg lik and prior from the modelInfo object
                now.logMargLik = nowInfo.logMargLik;
                now.logPrior = nowInfo.logPrior;
            }

            // decide acceptance:
            // for acceptance, the new model must be valid and the acceptance must be sampled
            if ((R_IsNaN(now.logMargLik) == FALSE) &&
                (unif_rand() <= exp(now.logMargLik - old.logMargLik + now.logPrior - old.logPrior + logPropRatio)))
            { // acceptance
                old = now;
            }
            else
            { // rejection
                now = old;
            }

            // so now definitely old == now, and we can
            // increment the associated sampling frequency.
            modelCache.incrementFrequency(now.modPar);

            // echo progress?
            if((++t % std::max(bookkeep.chainlength / 100, static_cast<PosLargeInt>(1)) == 0) &&
                bookkeep.verbose)
            {
                    Rprintf("-"); // display computation progress at each percent
            }
    }

    PutRNGstate(); // no RNs required anymore


    // normalize posterior probabilities and correct log marg lik and log prior
    const long double logNormConst = modelCache.getLogNormConstant();

    // get the nModels best models from the cache as an R list
    List ret = modelCache.getListOfBestModels(fpInfo,
                                              logNormConst,
                                              bookkeep);

    // set the attributes
    ret.attr("numVisited") = modelCache.size();
    ret.attr("inclusionProbs") = modelCache.getInclusionProbs(logNormConst, fpInfo.nFps, ucInfo.nUcGroups);
    ret.attr("logNormConst") = logNormConst;

    if (bookkeep.verbose){
        Rprintf("\nNumber of non-identifiable model proposals:     %d", bookkeep.nanCounter);
        Rprintf("\nNumber of total cached models:                  %d", modelCache.size());
        Rprintf("\nNumber of returned models:                      %d\n", Rf_length(ret));
    }

    // and return
    return ret;
}


// ***************************************************************************************************//


List
glmExhaustive(const DataValues& data,
              const FpInfo& fpInfo,
              const UcInfo& ucInfo,
              const FixInfo& fixInfo,
              Book& bookkeep,
              const GlmModelConfig& config,
              const GaussHermite& gaussHermite)
{
    // no map needed for exhaustive search, a set is the right thing:
    set<Model> orderedModels;

    // check that the set is large enough for all models
    if (orderedModels.max_size() < data.totalNumber)
        Rf_error("\nmodel space is too large - cannot compute every model\n");

    // start model
    ModelPar startModel(fpInfo.nFps);

    // bookkeeping

    // for computation of inclusion probs:
    // array of IndexSafeSum objects.
    //const int cgwp_length = fpInfo.nFps + ucInfo.nUcGroups;
    
    //IndexSafeSum cgwp[cgwp_length];
    //IndexSafeSum cgwp[fpInfo.nFps + ucInfo.nUcGroups];
    //bookkeep.covGroupWisePosteriors = cgwp;

    // now with vectors instead of arrays
    std::vector<IndexSafeSum> cgwp(fpInfo.nFps + ucInfo.nUcGroups);
    bookkeep.covGroupWisePosteriors = cgwp;
    
    
    // calculate the true null model if we have any fixed covariates,
    // otherwise it comes later
    if(fixInfo.nFixGroups != 0)
      computeGlm(startModel, orderedModels,
               data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);
    
    // add the fixed covariates to the model configuration
    IntSet s;
    for (unsigned int i = 0; i < fixInfo.nFixGroups; ++i) 
      s.insert(s.end(), i+1);
    startModel.fixPars = s;
      
    
    // start computation
    glmPermPars(0, startModel, orderedModels,
                data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);

    // we have finished.

    // now echo statistics:
    if (bookkeep.verbose)
    {
        Rprintf("\nActual number of possible models:  %d ",
                bookkeep.modelCounter);
        Rprintf("\nNumber of non-identifiable models: %d",
                bookkeep.nanCounter);
        Rprintf("\nNumber of saved possible models:   %d\n",
                orderedModels.size());
    }

    // normalize posterior probabilities of the models to return
    // (we do not know here the normalizing constant for the marginal likelihoods!)
    const long double logNormConst = bookkeep.modelLogPosteriors.logSumExp();

    // allocate the return list
    List ret(orderedModels.size());

    // and fill it:

    // first the single models
    R_len_t i = 0;
    for (set<Model>::const_reverse_iterator
            j = orderedModels.rbegin();
            j != orderedModels.rend();
            j++)
    {
        ret[i++] = j->convert2list(fpInfo,
                                   logNormConst,
                                   bookkeep);
    }

    // then some attributes:

    NumericVector inc(fpInfo.nFps + ucInfo.nUcGroups); // TODO should I add fix here
    for (R_len_t i = 0; i != inc.size(); ++i)
    {
        inc[i] = bookkeep.covGroupWisePosteriors[i].sumNormalizedExp(bookkeep.modelLogPosteriors, logNormConst);
    }
    ret.attr("inclusionProbs") = inc;
    ret.attr("numVisited") = static_cast<double>(bookkeep.modelCounter);
    ret.attr("logNormConst") = static_cast<double>(logNormConst);

    // finally return the list.
    return ret;
}

// ***************************************************************************************************//

// 21/11/2012: add tbf option
// 03/12/2012: add Cox regression with tbfs
// 03/07/2013: remove nullModelInfo and only get nullModelLogMargLik

// R call is:
//
//Ret <-
//    .External (cpp_glmBayesMfp,
//               data,
//               fpInfos,
//               ucInfos,
//               fixInfos,
//               searchConfig,
//               distribution,
//               options)
// [[Rcpp::export]]
SEXP
cpp_glmBayesMfp(List rcpp_data, List rcpp_fpInfos, List rcpp_ucInfos, List rcpp_fixInfos,
                List rcpp_searchConfig, List rcpp_distribution, List rcpp_options )
{
    // ----------------------------------------------------------------------------------
    // extract arguments
    // ----------------------------------------------------------------------------------

//     r_interface = CDR(r_interface);
//     List rcpp_data(CAR(r_interface));
// 
//     r_interface = CDR(r_interface);
//     List rcpp_fpInfos(CAR(r_interface));
// 
//     r_interface = CDR(r_interface);
//     List rcpp_ucInfos(CAR(r_interface));
// 
// 	  r_interface = CDR(r_interface);
// 	  List rcpp_fixInfos(CAR(r_interface));
// 
//     r_interface = CDR(r_interface);
//     List rcpp_searchConfig(CAR(r_interface));
// 
//     r_interface = CDR(r_interface);
//     List rcpp_distribution(CAR(r_interface));
// 
//     r_interface = CDR(r_interface);
//     List rcpp_options(CAR(r_interface));


    // ----------------------------------------------------------------------------------
    // unpack the R objects
    // ----------------------------------------------------------------------------------

    // data:
    const NumericMatrix n_x = rcpp_data["x"];
    const AMatrix x(n_x.begin(), n_x.nrow(),
                   n_x.ncol());

    const NumericMatrix n_xCentered = rcpp_data["xCentered"];
    const AMatrix xCentered(n_xCentered.begin(), n_xCentered.nrow(),
                           n_xCentered.ncol());

    const NumericVector n_y = rcpp_data["y"];
    const AVector y(n_y.begin(), n_y.size());

    const IntVector censInd = as<IntVector>(rcpp_data["censInd"]);

    // FP configuration:

    // vector of maximum fp degrees
    const PosIntVector fpmaxs = as<PosIntVector>(rcpp_fpInfos["fpmaxs"]);
    // corresponding vector of fp column indices
    const PosIntVector fppos = rcpp_fpInfos["fppos"];
    // corresponding vector of power set cardinalities
    const PosIntVector fpcards = rcpp_fpInfos["fpcards"];
    // names of fp terms
    const StrVector fpnames = rcpp_fpInfos["fpnames"];


    // UC configuration:

    const PosIntVector ucIndices = rcpp_ucInfos["ucIndices"];
    List rcpp_ucColList = rcpp_ucInfos["ucColList"];

    std::vector<PosIntVector> ucColList;
    for (R_len_t i = 0; i != rcpp_ucColList.length(); ++i)
    {
        ucColList.push_back(as<PosIntVector>(rcpp_ucColList[i]));
    }


	// Fixed Covariate configuration:

	const PosIntVector fixIndices = rcpp_fixInfos["fixIndices"];
	List rcpp_fixColList = rcpp_fixInfos["fixColList"];

	std::vector<PosIntVector> fixColList;
	for (R_len_t i = 0; i != rcpp_fixColList.length(); ++i)
	{
		fixColList.push_back(as<PosIntVector>(rcpp_fixColList[i]));
	}


    
    // model search configuration:

    const double totalNumber = as<double>(rcpp_searchConfig["totalNumber"]);
    const PosInt nModels = as<PosInt>(rcpp_searchConfig["nModels"]);
    const bool empiricalBayes = as<bool>(rcpp_searchConfig["empiricalBayes"]);
    const bool useFixedg = as<bool>(rcpp_searchConfig["useFixedg"]);
    const bool doSampling = as<bool>(rcpp_searchConfig["doSampling"]);
    const double chainlength = as<double>(rcpp_searchConfig["chainlength"]);
    const PosInt nCache = as<PosInt>(rcpp_searchConfig["nCache"]);
    const double largeVariance = as<double>(rcpp_searchConfig["largeVariance"]);
    const bool useBfgs = as<bool>(rcpp_searchConfig["useBfgs"]);
    const bool useFixedc = as<bool>(rcpp_searchConfig["useFixedc"]);

    // there might be a single model configuration saved in the searchConfig:
    bool onlyComputeModelsInList;
    try {
        onlyComputeModelsInList = true;
        List rcpp_modelConfigs = rcpp_searchConfig["modelConfigs"];
    } catch (std::exception& e) {
        onlyComputeModelsInList = false;
    }
    // we need to try-catch it because there might be no element "modelConfigs" in
    // rcpp_searchConfig in which case Rcpp throws an exception.

    // distributions info:

    const std::string modelPrior = as<std::string>(rcpp_distribution["modelPrior"]);
    const bool doGlm = as<bool>(rcpp_distribution["doGlm"]);
    const bool tbf = as<bool>(rcpp_distribution["tbf"]);
    const double nullModelLogMargLik = as<double>(rcpp_distribution["nullModelLogMargLik"]);
    const double nullModelDeviance = as<double>(rcpp_distribution["nullModelDeviance"]);
    const double fixedg = as<double>(rcpp_distribution["fixedg"]);
    const double empiricalMean = as<double>(rcpp_distribution["yMean"]);
    const bool empiricalgPrior = as<bool>(rcpp_distribution["empiricalgPrior"]);
    
    S4 rcpp_gPrior = rcpp_distribution["gPrior"];
    List rcpp_family = rcpp_distribution["family"];

    // other options:

    const bool verbose = as<bool>(rcpp_options["verbose"]);
    const bool debug = as<bool>(rcpp_options["debug"]);
#ifdef _OPENMP
    const bool useOpenMP = as<bool>(rcpp_options["useOpenMP"]);
#endif
    const GaussHermite gaussHermite(as<List>(rcpp_options["gaussHermite"]));
    const bool higherOrderCorrection = as<bool>(rcpp_options["higherOrderCorrection"]);


    // ----------------------------------------------------------------------------------
    // further process input information
    // ----------------------------------------------------------------------------------

    // data:

    // only the intercept is always included, that is fixed, in the model
    IntSet fixedCols;
    fixedCols.insert(1); //vestigial code that isn't part of the new support for fixed vars

    const DataValues data(x, xCentered, y, censInd, totalNumber, fixedCols);

    // FP configuration:
    const FpInfo fpInfo(fpcards, fppos, fpmaxs, fpnames, x);

    // UC configuration:

    // determine sizes of the UC groups, and the total size == maximum size reached together by all
    // UC groups.
    PosIntVector ucSizes;
    PosInt maxUcDim = 0;
    for (vector<PosIntVector>::const_iterator cols = ucColList.begin(); cols != ucColList.end(); ++cols)
    {
        PosInt thisSize = cols->size();

        maxUcDim += thisSize;
        ucSizes.push_back(thisSize);
    }
    const UcInfo ucInfo(ucSizes, maxUcDim, ucIndices, ucColList);


	// Fix configuration:

	// determine sizes of the fixed covariate groups, and the total size == maximum size reached together by all
	// fix groups.
	PosIntVector fixSizes;
	PosInt maxFixDim = 0;
	for (vector<PosIntVector>::const_iterator cols = fixColList.begin(); cols != fixColList.end(); ++cols)
	{
		PosInt thisSize = cols->size();

		maxFixDim += thisSize;
		fixSizes.push_back(thisSize);
	}
	const FixInfo fixInfo(fixSizes, maxFixDim, fixIndices, fixColList);




    // model search configuration:
    Book bookkeep(tbf,
                  doGlm,
                  empiricalBayes,
                  useFixedg,
                  useFixedc,
                  chainlength,
                  doSampling,
                  verbose,
                  modelPrior,
                  nModels,
                  nCache,
                  largeVariance,
                  useBfgs,
                  debug,
                  higherOrderCorrection);

    // model configuration:
    const GlmModelConfig config(rcpp_family, nullModelLogMargLik, nullModelDeviance, fixedg, rcpp_gPrior,
                                data.response, bookkeep.debug, bookkeep.useFixedc, empiricalMean,
                                empiricalgPrior);

    // use only one thread if we do not want to use openMP.
#ifdef _OPENMP
    if(! useOpenMP)
    {
        omp_set_num_threads(1);
    } else {
        // else use all available cpu's.
        omp_set_num_threads(omp_get_num_procs());
    }
#endif

    // ----------------------------------------------------------------------------------
    // now either compute only one Model, do model sampling or do an exhaustive search
    // ----------------------------------------------------------------------------------

    if(onlyComputeModelsInList)
    {
        return glmModelsInList(data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite, as<List>(rcpp_searchConfig["modelConfigs"]));
    }
    else if(doSampling)
    {
        return glmSampling(data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);
    }
    else
    {
        return glmExhaustive(data, fpInfo, ucInfo, fixInfo, bookkeep, config, gaussHermite);
    }
}

// ***************************************************************************************************//

// End of file.
