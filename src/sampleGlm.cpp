/*
 * sampleGlm.cpp
 *
 *  Created on: 10.12.2009
 *      Author: daniel
 *      
 *      13/07/2015 Replace assert() with Rccp:Stop()
 */

#include <rcppExport.h>
#include <combinatorics.h>
#include <dataStructure.h>
#include <types.h>
#include <iwls.h>
#include <design.h>
#include <coxfit.h>
#include <bfgs.h>
#include <optimize.h>
#include <fpUcHandling.h>
#include <linalgInterface.h>
//#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ***************************************************************************************************//

struct MarginalZ
{
    MarginalZ(const RFunction& logDens,
              const RFunction& gen) :
                  logDens(logDens),
                  gen(gen)
                  {
                  }


    const RFunction logDens;
    const RFunction gen;
};

// ***************************************************************************************************//

struct Options
{
    Options(bool estimateMargLik,
            bool verbose,
            bool debug,
            bool isNullModel,
            bool useFixedZ,
            bool tbf,
            bool doGlm,
            PosInt iterations,
            PosInt burnin,
            PosInt step) :
        estimateMargLik(estimateMargLik), verbose(verbose),
                debug(debug), isNullModel(isNullModel), useFixedZ(useFixedZ),
                tbf(tbf), doGlm(doGlm),
                nSamples(ceil((iterations - burnin) * 1.0 / step)),
                iterations(iterations), burnin(burnin),
                step(step)
    {
    }

    const bool estimateMargLik;
    const bool verbose;
    const bool debug;
    const bool isNullModel;
    const bool useFixedZ;
    const bool tbf;
    const bool doGlm;

    const PosInt nSamples;
    const PosInt iterations;
    const PosInt burnin;
    const PosInt step;
};

// ***************************************************************************************************//

// this is just for being sure that the memory is returned correctly even in case
// of an error or user interrupt from R.
struct Fitter
{
    Fitter() :
        iwlsObject(0),
        coxfitObject(0)
    {}

    ~Fitter()
    {
        delete iwlsObject;
        delete coxfitObject;
    }

    Iwls * iwlsObject;
    Coxfit * coxfitObject;
};



// ***************************************************************************************************//


class Mcmc
{

public:
    // ctr
    Mcmc(const MarginalZ& marginalz, PosInt nObs, PosInt nCoefs) :
        sample(nCoefs),
        proposalInfo(nObs, nCoefs),
        marginalz(marginalz)
        {
        }

    // the current parameter sample
    Parameter sample;

    // the unnormalized log posterior of this sample
    double logUnPosterior;

    // info about the normal proposal distribution given the sampled z
    IwlsResults proposalInfo;

    // compute the log of the normalized proposal density when the z log density is provided
    // normalize correctly in order not to get problems (perhaps) with the Chib-Jeliazkov estimate
    // computation.
    double
    computeLogProposalDens() const
    {
        // Be careful: qFactor is in fact lower-triangular, so a simple multiplication would fail!
        // use instead directly a BLAS routine for this multiplication.
        AVector tmp = sample.coefs - proposalInfo.coefs;
        trmv(false, true, proposalInfo.qFactor, tmp);

        return 0.5 * (proposalInfo.logPrecisionDeterminant - arma::dot(tmp, tmp)) -
               M_LN_SQRT_2PI * proposalInfo.qFactor.n_rows +
               marginalz.logDens(sample.z);
    }

    // non-default assignment operator
    Mcmc&
    operator=(const Mcmc& rhs)
    {
        if(this == &rhs)
        {
            return *this;
        }
        else
        {
            sample = rhs.sample;
            logUnPosterior = rhs.logUnPosterior;
            proposalInfo = rhs.proposalInfo;

            return *this;
        }
    }


private:
    // the marginal z info: same for all Mcmc objects,
    // therefore it is not assigned by the assignment operator
    const MarginalZ marginalz;
    // important: copy the object, because otherwise (reference/pointer)
    // we are not sure that the functions are still available if we do not use "new"
};

// ***************************************************************************************************//


struct Samples
{
    // constructor: basically allocates beta matrix
    Samples(PosInt nCoefs, PosInt nSamples) :
        coefsSamples(nCoefs, nSamples),
        nSaved(0)
        {
        }

    // save a sample consisting of coefs and z
    void
    storeParameters(const Parameter& sample)
    {
        coefsSamples.col(nSaved++) = sample.coefs;
        zSamples.push_back(sample.z);
    }

    // save terms for marginal likelihood estimate
    void
    storeMargLikTerms(double num, double denom)
    {
        numerator.push_back(num);
        denominator.push_back(denom);
    }

    // output everything to an R list
    List
    convert2list() const;

private:
    // nCoefs x nSamples:
    AMatrix coefsSamples;

    // counts the number of saved parameters, so that we know where to store the
    // next coefficients vector
    PosInt nSaved;

    // is gradually extended:
    MyDoubleVector zSamples;

    // possibly stays empty if not required by the user:
    // the numerator and denominator terms for the marginal likelihood estimate
    MyDoubleVector numerator;
    MyDoubleVector denominator;
};

List
Samples::convert2list() const
{
    return List::create(_["coefficients"] = coefsSamples,
                        _["z"] = zSamples,
                        _["margLikNumerator"] = numerator,
                        _["margLikDenominator"] = denominator);
}


// ***************************************************************************************************//


// get a vector with normal variates from N(mean, sd^2)
AVector
drawNormalVariates(PosInt n, double mean, double sd)
{
    AVector ret(n);

    // use R's random number generator
    GetRNGstate();

    for (PosInt i = 0; i < n; ++i)
    {
        ret(i) = Rf_rnorm(mean, sd);
    }

    // no RNs required anymore
    PutRNGstate();

    return ret;
}

// draw a single random normal vector from N(mean, (precisionCholeskyFactor * t(precisionCholeskyFactor))^(-1))
AVector
drawNormalVector(const AVector& mean,
                 const AMatrix& precisionCholeskyFactor)
{
    // get vector from N(0, I)
    AVector w = drawNormalVariates(mean.n_rows, // as many normal variates as required by the dimension.
                                   0.0,
                                   1.0);

    // then solve L' * ret = w, and overwrite w with the result:
    trs(false,
        true,
        precisionCholeskyFactor,
        w);

    // return the shifted vector
    return (w + mean);
}


// draw a single uniform random variable:
// be careful with the seed because the z generator function also uses it (via R)
double
unif()
{
    GetRNGstate();

    double ret = unif_rand();

    PutRNGstate();

    return ret;
}


// ***************************************************************************************************//


// R call is:
//
//    samples <- .External(cpp_sampleGlm,
//                         model,
//                         attrs$data,
//                         attrs$fpInfos,
//                         attrs$ucInfos,
//                         attrs$fixInfos,
//                         attrs$distribution,
//                         newdata,
//                         options,
//                         marginalz)




// [[Rcpp::export]]
SEXP
cpp_sampleGlm(      List rcpp_model,List rcpp_data, List rcpp_fpInfos, List rcpp_ucInfos,
                    List rcpp_fixInfos, List rcpp_distribution, List rcpp_searchConfig,
                    List rcpp_options, List rcpp_marginalz)
{
    // ----------------------------------------------------------------------------------
    // extract arguments
    // ----------------------------------------------------------------------------------

    // r_interface = CDR(r_interface);
    // List rcpp_model(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_data(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_fpInfos(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_ucInfos(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_fixInfos(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_distribution(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_searchConfig(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_options(CAR(r_interface));
    // 
    // r_interface = CDR(r_interface);
    // List rcpp_marginalz(CAR(r_interface));

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
    
    // fixed covariate configuration:
    
    const PosIntVector fixIndices = rcpp_fixInfos["fixIndices"];
    List rcpp_fixColList = rcpp_fixInfos["fixColList"];
    
    std::vector<PosIntVector> fixColList;
    for (R_len_t i = 0; i != rcpp_fixColList.length(); ++i)
    {
      fixColList.push_back(as<PosIntVector>(rcpp_fixColList[i]));
    }
    
    
    

    // distributions info:

    const double nullModelLogMargLik = as<double>(rcpp_distribution["nullModelLogMargLik"]);
    const double nullModelDeviance = as<double>(rcpp_distribution["nullModelDeviance"]);
    S4 rcpp_gPrior = rcpp_distribution["gPrior"];
    List rcpp_family = rcpp_distribution["family"];
    const bool tbf = as<bool>(rcpp_distribution["tbf"]);
    const bool doGlm = as<bool>(rcpp_distribution["doGlm"]);
    
    const double empiricalMean = as<double>(rcpp_distribution["yMean"]);
    const bool empiricalgPrior = as<bool>(rcpp_distribution["empiricalgPrior"]);
    
    // model search configuration:
    const bool useFixedc = as<bool>(rcpp_searchConfig["useFixedc"]);
    
    // options:

    const bool estimateMargLik = as<bool>(rcpp_options["estimateMargLik"]);
    const bool verbose = as<bool>(rcpp_options["verbose"]);
    const bool debug = as<bool>(rcpp_options["debug"]);
    const bool isNullModel = as<bool>(rcpp_options["isNullModel"]);
    const bool useFixedZ = as<bool>(rcpp_options["useFixedZ"]);
    const double fixedZ = as<double>(rcpp_options["fixedZ"]);
#ifdef _OPENMP
    const bool useOpenMP = as<bool>(rcpp_options["useOpenMP"]);
#endif

    S4 rcpp_mcmc = rcpp_options["mcmc"];
    const PosInt iterations = rcpp_mcmc.slot("iterations");
    const PosInt burnin = rcpp_mcmc.slot("burnin");
    const PosInt step = rcpp_mcmc.slot("step");


    // z density stuff:

    const RFunction logMarginalZdens(as<SEXP>(rcpp_marginalz["logDens"]));
    const RFunction marginalZgen(as<SEXP>(rcpp_marginalz["gen"]));


    // ----------------------------------------------------------------------------------
    // further process arguments
    // ----------------------------------------------------------------------------------

    // data:

     // only the intercept is always included, that is fixed, in the model
     IntSet fixedCols;
     fixedCols.insert(1);

     // totalnumber is set to 0 because we do not care about it.
     const DataValues data(x, xCentered, y, censInd, 0, fixedCols);

     // FP configuration:
     const FpInfo fpInfo(fpcards, fppos, fpmaxs, fpnames, x);

     // UC configuration:

     // determine sizes of the UC groups, and the total size == maximum size reached together by all
     // UC groups.
     PosIntVector ucSizes;
     PosInt maxUcDim = 0;
     for (std::vector<PosIntVector>::const_iterator cols = ucColList.begin(); cols != ucColList.end(); ++cols)
     {
         PosInt thisSize = cols->size();

         maxUcDim += thisSize;
         ucSizes.push_back(thisSize);
     }
     const UcInfo ucInfo(ucSizes, maxUcDim, ucIndices, ucColList);

  
     // fix configuration:  
     
     // determine sizes of the fix groups, and the total size == maximum size reached together by all
     // UC groups.
     PosIntVector fixSizes;
     PosInt maxFixDim = 0;
     for (std::vector<PosIntVector>::const_iterator cols = fixColList.begin(); cols != fixColList.end(); ++cols)
     {
       PosInt thisSize = cols->size();
       
       maxFixDim += thisSize;
       fixSizes.push_back(thisSize);
     }
     const FixInfo fixInfo(fixSizes, maxFixDim, fixIndices, fixColList);
     
     
     
     
     
     // model configuration:
     GlmModelConfig config(rcpp_family, nullModelLogMargLik, nullModelDeviance, exp(fixedZ), rcpp_gPrior,
                           data.response, debug, useFixedc, empiricalMean, empiricalgPrior);


     // model config/info:
     const Model thisModel(ModelPar(rcpp_model["configuration"],
                                   fpInfo),
                          GlmModelInfo(as<List>(rcpp_model["information"])));


     // the options
     const Options options(estimateMargLik,
                           verbose,
                           debug,
                           isNullModel,
                           useFixedZ,
                           tbf,
                           doGlm,
                           iterations,
                           burnin,
                           step);

     // marginal z stuff
     const MarginalZ marginalZ(logMarginalZdens,
                               marginalZgen);


     // use only one thread if we do not want to use openMP.
#ifdef _OPENMP
     if(! useOpenMP)
     {
         omp_set_num_threads(1);
     } else {
         omp_set_num_threads(omp_get_num_procs());
     }
#endif


     // ----------------------------------------------------------------------------------
     // prepare the sampling
     // ----------------------------------------------------------------------------------

     Fitter fitter;
     int nCoefs;

     if(options.doGlm)
     {
         // construct IWLS object, which can be used for all IWLS stuff,
         // and also contains the design matrix etc
         fitter.iwlsObject = new Iwls(thisModel.par,
                                      data,
                                      fpInfo,
                                      ucInfo,
                                      fixInfo,
                                      config,
                                      config.linPredStart,
                                      options.useFixedZ,
                                      EPS,
                                      options.debug,
                                      options.tbf);

         nCoefs = fitter.iwlsObject->nCoefs;

         // check that we have the same answer about the null model as R
         //assert(fitter.iwlsObject->isNullModel == options.isNullModel);
         if(fitter.iwlsObject->isNullModel != options.isNullModel){
           Rcpp::stop("sampleGlm.cpp:cpp_sampleGlm: isNullModel != options.isNullModel");
         } 
     }
     else
     {
         AMatrix design = getDesignMatrix(thisModel.par, data, fpInfo, ucInfo, fixInfo, false);
         fitter.coxfitObject = new Coxfit(data.response,
                                          data.censInd,
                                          design,
                                          config.weights,
                                          config.offsets,
                                          1);

         // the number of coefficients (here it does not include the intercept!!)
         nCoefs = design.n_cols;

         // check that we do not have a null model here:
         // assert(nCoefs > 0);
         if(nCoefs <= 0){
           Rcpp::stop("sampleGlm.cpp:cpp_sampleGlm: nCoefs <= 0");
         } 
     }


     // allocate sample container
     Samples samples(nCoefs, options.nSamples);

     // count how many proposals we have accepted:
     PosInt nAccepted(0);

     // at what z do we start?
     double startZ = useFixedZ ? fixedZ : thisModel.info.zMode;

     // start container with current things
     Mcmc now(marginalZ, data.nObs, nCoefs);

     if(doGlm)
     {
         // get the mode for beta given the mode of the approximated marginal posterior as z
         // if TBF approach is used, this will be the only time the IWLS is used,
         // because we only need the MLE and the Cholesky factor of its
         // precision matrix estimate, which do not depend on z.
         PosInt iwlsIterations = fitter.iwlsObject->startWithNewLinPred(40,
                                                                        // this is the corresponding g
                                                                        exp(startZ),
                                                                        // and the start value for the linear predictor is taken from the Glm model config
                                                                        config.linPredStart);

         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Initial IWLS for high density point finished after %d iterations",
                     iwlsIterations);
         }

         // this is the current proposal info:
         now.proposalInfo = fitter.iwlsObject->getResults();

         // and this is the current parameters sample:
         now.sample = Parameter(now.proposalInfo.coefs,
                                startZ);

         if(options.tbf)
         {
             // we will not compute this in the TBF case:
             now.logUnPosterior = R_NaReal;

             // start to compute the variance of the intercept parameter:

             // here the inverse cholesky factor of the precision matrix will
             // be stored. First, it's the identity matrix.
             AMatrix inverseQfactor = arma::eye(now.proposalInfo.qFactor.n_rows,
                                                now.proposalInfo.qFactor.n_cols);

             // do the inversion
             trs(false,
                 false,
                 now.proposalInfo.qFactor,
                 inverseQfactor);

             // now we can compute the variance of the intercept estimate:
             const AVector firstCol = inverseQfactor.col(0);
             const double interceptVar = arma::dot(firstCol, firstCol);

             // ok, now alter the qFactor appropriately to reflect the
             // independence assumption between the intercept estimate
             // and the other coefficients estimates
             now.proposalInfo.qFactor.col(0) = arma::zeros<AVector>(now.proposalInfo.qFactor.n_rows);
             now.proposalInfo.qFactor(0, 0) = sqrt(1.0 / interceptVar);
         }
         else
         {
             // compute the (unnormalized) log posterior of the proposal
             now.logUnPosterior = fitter.iwlsObject->computeLogUnPosteriorDens(now.sample);
         }
     }
     else
     {
         PosInt coxfitIterations = fitter.coxfitObject->fit();
         CoxfitResults coxResults = fitter.coxfitObject->finalizeAndGetResults();
         fitter.coxfitObject->checkResults();

         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Cox fit finished after %d iterations",
                     coxfitIterations);
         }

         // we will not compute this in the TBF case:
         now.logUnPosterior = R_NaReal;

         // compute the Cholesky factorization of the covariance matrix
         int info = potrf(false,
                          coxResults.imat);

         // check that all went well
         if(info != 0)
         {
             std::ostringstream stream;
             stream << "dpotrf(coxResults.imat) got error code " << info << "in sampleGlm";
             throw std::domain_error(stream.str().c_str());
         }

         // compute the precision matrix, using the Cholesky factorization
         // of the covariance matrix
         now.proposalInfo.qFactor = arma::eye(now.proposalInfo.qFactor.n_rows,
                                              now.proposalInfo.qFactor.n_cols);
         info = potrs(false,
                      coxResults.imat,
                      now.proposalInfo.qFactor);

         // check that all went well
         if(info != 0)
         {
             std::ostringstream stream;
             stream << "dpotrs(coxResults.imat, now.proposalInfo.qFactor) got error code " << info << "in sampleGlm";
             throw std::domain_error(stream.str().c_str());
         }

         // compute the Cholesky factorization of the precision matrix
         info = potrf(false,
                      now.proposalInfo.qFactor);

         // check that all went well
         if(info != 0)
         {
             std::ostringstream stream;
             stream << "dpotrf(now.proposalInfo.qFactor) got error code " << info << "in sampleGlm";
             throw std::domain_error(stream.str().c_str());
         }

         // the MLE of the coefficients
         now.proposalInfo.coefs = coxResults.coefs;
     }

     // so the parameter object "now" is then also the high density point
     // required for the marginal likelihood estimate:
     const Mcmc highDensityPoint(now);

     // we accept this starting value, so initialize "old" with the same ones
     Mcmc old(now);

     // ----------------------------------------------------------------------------------
     // start sampling
     // ----------------------------------------------------------------------------------

     // echo debug-level message?
     if(options.debug)
     {
         if(tbf)
         {
             Rprintf("\ncpp_sampleGlm: Starting MC simulation");
         }
         else
         {
             Rprintf("\ncpp_sampleGlm: Starting MCMC loop");
         }
     }


     // i_iter starts at 1 !!
     for(PosInt i_iter = 1; i_iter <= options.iterations; ++i_iter)
     {
         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Starting iteration no. %d", i_iter);
         }

         // ----------------------------------------------------------------------------------
         // store the proposal
         // ----------------------------------------------------------------------------------

         // sample one new log covariance factor z (other arguments than 1 are not useful
         // with the current setup of the RFunction wrapper class)
         now.sample.z = marginalZ.gen(1);

         if(options.tbf)
         {
             if(options.isNullModel)
             {
                 // note that we do not encounter this in the Cox case
                 // assert(options.doGlm);
                 if(!options.doGlm){
                   Rcpp::stop("sampleGlm.cpp:cpp_sampleGlm: options.doGlm should be TRUE");
                 } 

                 // draw the proposal coefs, which is here just the intercept
                 now.sample.coefs = drawNormalVector(now.proposalInfo.coefs,
                                                     now.proposalInfo.qFactor);

             }
             else
             {   // here we have at least one non-intercept coefficient

                 // get vector from N(0, I)
                 AVector w = drawNormalVariates(now.proposalInfo.coefs.n_elem,
                                                0.0,
                                                1.0);

                 // then solve L' * ret = w, and overwrite w with the result:
                 trs(false,
                     true,
                     now.proposalInfo.qFactor,
                     w);

                 // compute the shrinkage factor t = g / (g + 1)
                 const double g = exp(now.sample.z);
                
                //Previously used g directly, but if g=inf we need to use the limit
                 // const double shrinkFactor = g / (g + 1.0);
                const double shrinkFactor = std::isinf(g) ? 1 : g / (g + 1.0);
                 
                 // scale the variance of the non-intercept coefficients
                 // with this factor.
                 // In the Cox case: no intercept present, so scale everything
                 int startCoef = options.doGlm ? 1 : 0;

                 w.rows(startCoef, w.n_rows - 1) *= sqrt(shrinkFactor);

                 // also scale the mean of the non-intercept coefficients
                 // appropriately:
                 // In the Cox case: no intercept present, so scale everything
                 now.sample.coefs = now.proposalInfo.coefs;
                 now.sample.coefs.rows(startCoef, now.sample.coefs.n_rows - 1) *= shrinkFactor;

                 // so altogether we have:
                 now.sample.coefs += w;
             }
             ++nAccepted;
         }
         else // the generalized hyper-g prior case
         {
             // do 1 IWLS step, starting from the last linear predictor and the new z
             // (here the return value is not very interesting, as it must be 1)
             fitter.iwlsObject->startWithNewCoefs(1,
                                                  exp(now.sample.z),
                                                  now.sample.coefs);

             // get the results
             now.proposalInfo = fitter.iwlsObject->getResults();

             // draw the proposal coefs:
             now.sample.coefs = drawNormalVector(now.proposalInfo.coefs,
                                                 now.proposalInfo.qFactor);

             // compute the (unnormalized) log posterior of the proposal
             now.logUnPosterior = fitter.iwlsObject->computeLogUnPosteriorDens(now.sample);

             // ----------------------------------------------------------------------------------
             // get the reverse jump normal density
             // ----------------------------------------------------------------------------------

             // copy the old Mcmc object
             Mcmc reverse(old);

             // do again 1 IWLS step, starting from the sampled linear predictor and the old z
             fitter.iwlsObject->startWithNewCoefs(1,
                                          exp(reverse.sample.z),
                                          now.sample.coefs);

             // get the results for the reverse jump Gaussian:
             // only the proposal has changed in contrast to the old container,
             // the sample stays the same!
             reverse.proposalInfo = fitter.iwlsObject->getResults();


             // ----------------------------------------------------------------------------------
             // compute the proposal density ratio
             // ----------------------------------------------------------------------------------

             // first the log of the numerator, i.e. log(f(old | new)):
             double logProposalRatioNumerator = reverse.computeLogProposalDens();

             // second the log of the denominator, i.e. log(f(new | old)):
             double logProposalRatioDenominator = now.computeLogProposalDens();

             // so the log proposal density ratio is
             double logProposalRatio = logProposalRatioNumerator - logProposalRatioDenominator;

             // ----------------------------------------------------------------------------------
             // compute the posterior density ratio
             // ----------------------------------------------------------------------------------

             double logPosteriorRatio = now.logUnPosterior - old.logUnPosterior;

             // ----------------------------------------------------------------------------------
             // accept or reject proposal
             // ----------------------------------------------------------------------------------

             double acceptanceProb = exp(logPosteriorRatio + logProposalRatio);

             if(unif() < acceptanceProb)
             {
                 old = now;

                 ++nAccepted;
             }
             else
             {
                 now = old;
             }
         }

         // ----------------------------------------------------------------------------------
         // store the sample?
         // ----------------------------------------------------------------------------------

         // if the burnin was passed and we are at a multiple of step beyond that, then store
         // the sample.
         if((i_iter > options.burnin) &&
            (((i_iter - options.burnin) % options.step) == 0))
         {
             // echo debug-level message
             if(options.debug)
             {
                 Rprintf("\ncpp_sampleGlm: Storing samples of iteration no. %d", i_iter);
             }

             // store the current parameter sample
             samples.storeParameters(now.sample);

             // ----------------------------------------------------------------------------------
             // compute marginal likelihood terms
             // ----------------------------------------------------------------------------------

             // compute marginal likelihood terms and save them?
             // (Note that the tbf bool is just for safety here,
             // the R function sampleGlm will set estimateMargLik to FALSE
             // when tbf is TRUE.)
             if(options.estimateMargLik && (! options.tbf))
             {
                 // echo debug-level message?
                 if(options.debug)
                 {
                     Rprintf("\ncpp_sampleGlm: Compute marginal likelihood estimation terms");
                 }

                 // ----------------------------------------------------------------------------------
                 // compute next term for the denominator
                 // ----------------------------------------------------------------------------------

                 // draw from the high density point proposal distribution
                 Mcmc denominator(highDensityPoint);
                 denominator.sample.z = marginalZ.gen(1);

                 fitter.iwlsObject->startWithNewLinPred(1,
                                                exp(denominator.sample.z),
                                                highDensityPoint.proposalInfo.linPred);

                 denominator.proposalInfo = fitter.iwlsObject->getResults();

                 denominator.sample.coefs = drawNormalVector(denominator.proposalInfo.coefs,
                                                             denominator.proposalInfo.qFactor);

                 // get posterior density of the sample
                 denominator.logUnPosterior = fitter.iwlsObject->computeLogUnPosteriorDens(denominator.sample);

                 // get the proposal density at the sample
                 double denominator_logProposalDensity = denominator.computeLogProposalDens();

                 // then the reverse stuff:
                 // first we copy again the high density point
                 Mcmc revDenom(highDensityPoint);

                 // but choose the new sampled coefficients as starting point
                 fitter.iwlsObject->startWithNewCoefs(1,
                                              exp(revDenom.sample.z),
                                              denominator.sample.coefs);
                 revDenom.proposalInfo = fitter.iwlsObject->getResults();

                 // so the reverse proposal density is
                 double revDenom_logProposalDensity = revDenom.computeLogProposalDens();


                 // so altogether the next term for the denominator is the following acceptance probability
                 double denominatorTerm = denominator.logUnPosterior - highDensityPoint.logUnPosterior +
                                          revDenom_logProposalDensity - denominator_logProposalDensity;
                 denominatorTerm = exp(fmin(0.0, denominatorTerm));

                 // ----------------------------------------------------------------------------------
                 // compute next term for the numerator
                 // ----------------------------------------------------------------------------------

                 // compute the proposal density of the current sample starting from the high density point
                 Mcmc numerator(now);

                 fitter.iwlsObject->startWithNewLinPred(1,
                                                exp(numerator.sample.z),
                                                highDensityPoint.proposalInfo.linPred);
                 numerator.proposalInfo = fitter.iwlsObject->getResults();

                 double numerator_logProposalDensity = numerator.computeLogProposalDens();

                 // then compute the reverse proposal density of the high density point when we start from the current
                 // sample
                 Mcmc revNum(highDensityPoint);

                 fitter.iwlsObject->startWithNewCoefs(1,
                                              exp(revNum.sample.z),
                                              now.sample.coefs);
                 revNum.proposalInfo = fitter.iwlsObject->getResults();

                 double revNum_logProposalDensity = revNum.computeLogProposalDens();

                 // so altogether the next term for the numerator is the following guy:
                 double numeratorTerm = exp(fmin(revNum_logProposalDensity,
                                                 highDensityPoint.logUnPosterior - now.logUnPosterior +
                                                 numerator_logProposalDensity));

                 // ----------------------------------------------------------------------------------
                 // finally store both terms
                 // ----------------------------------------------------------------------------------

                 samples.storeMargLikTerms(numeratorTerm, denominatorTerm);

             }
         }

         // ----------------------------------------------------------------------------------
         // echo progress?
         // ----------------------------------------------------------------------------------

         // echo debug-level message?
         if(options.debug)
         {
             Rprintf("\ncpp_sampleGlm: Finished iteration no. %d", i_iter);
         }

         if((i_iter % std::max(static_cast<int>(options.iterations / 100), 1) == 0) &&
             options.verbose)
         {
             // display computation progress at each percent
             Rprintf("-");

         } // end echo progress

     } // end MCMC loop


     // echo debug-level message?
     if(options.debug)
     {
         if(tbf)
         {
             Rprintf("\ncpp_sampleGlm: Finished MC simulation");
         }
         else
         {
             Rprintf("\ncpp_sampleGlm: Finished MCMC loop");
         }
     }


     // ----------------------------------------------------------------------------------
     // build up return list for R and return that.
     // ----------------------------------------------------------------------------------

     return List::create(_["samples"] = samples.convert2list(),
                         _["nAccepted"] = nAccepted,
                         _["highDensityPointLogUnPosterior"] = highDensityPoint.logUnPosterior);

} // end cpp_sampleGlm

// ***************************************************************************************************//

// End of sampleGlm.cpp

