#include <set>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include <dataStructure.h>
#include <sum.h>
#include <functionWraps.h>
#include <fpUcHandling.h>

#include <rcppExport.h>

using std::lexicographical_compare;
using std::accumulate;
using std::exp;
using std::log;

using namespace Rcpp;

// ***************************************************************************************************//

// safeSum //

void SafeSum::add(const long double &val)
{
        vals.push_back(val);
}

long double SafeSum::sum()
{
        long double ret = modified_deflation(vals);
        return ret;
}

// compute the log of the sum of the exp of the elements using accurate algorithm,
// and avoiding infinite contributions.
long double SafeSum::logSumExp()
{
    // be sure that there is at least 1 value in "vals"
    if(vals.empty())
    {
        return R_NaN;
    }

    // the maximum of the log contributions is:
    long double maxLogContrib = *std::max_element(vals.begin(), vals.end());

    // now compute the constant which is added to all log contributions,
    // in order to avoid infinite contributions and at the same time use
    // the whole number space (i.e. possibly avoid zero contributions)
    long double constant = log(LDBL_MAX) - 100.0L - maxLogContrib;
    // 100 is for safety.

    // so now the contributions, offset by the constant
    LongDoubleVector expVals;
    for(LongDoubleVector::const_iterator
            l = vals.begin();
            l != vals.end();
            ++l)
    {
        expVals.push_back(exp(*l + constant));
    }

    // the result is the log of the sum, corrected with the constant:
    long double ret = log(modified_deflation(expVals)) - constant;
    return ret;
}

long double SafeSum::simpleSum()
{
        long double ret = 0.0;
        for(LongDoubleVector::const_iterator
                v = vals.begin();
                v != vals.end();
                ++v)
        {
            ret += *v;
        }
        return ret;
}

// ***************************************************************************************************//

// indexSafeSum //

void IndexSafeSum::add(const Ind& ind)
{
        indices.insert(ind);
}

long double
IndexSafeSum::sumNormalizedExp(const SafeSum& s, long double logNormConst) const
{
    LongDoubleVector tempVec;
    for (IndSet::const_iterator i = indices.begin(); i != indices.end(); i++)
    {
        tempVec.push_back(exp(s.vals.at(* i) - logNormConst));
    }
    return modified_deflation(tempVec);
}

// ***************************************************************************************************//


// Book //

// 21/11/2012: add tbf option

// constructor which checks the chainlength
Book::Book(bool tbf,
           bool doGlm,
           bool empiricalBayes,
           bool useFixedg,
           bool useFixedc,
           double cl,
           bool doSampling,
           bool verbose,
           std::string modelPrior,
           PosInt nModels,
           PosInt nCache,
           double largeVariance,
           bool useBfgs,
           bool debug,
           bool higherOrderCorrection) :
           modelCounter(0),
                nanCounter(0),
                tbf(tbf),
                doGlm(doGlm),
                empiricalBayes(empiricalBayes),
                useFixedg(useFixedg),
                useFixedc(useFixedc),
                doSampling(doSampling),
                verbose(verbose),
                modelPrior(modelPrior),
                nModels(nModels),
                nCache(nCache),
                largeVariance(largeVariance),
                useBfgs(useBfgs),
                debug(debug),
                higherOrderCorrection(higherOrderCorrection)
{
    if (doSampling)
    {
        // check if chainlength is too large
        if (ULLONG_MAX < cl)
        {
            Rf_warning("\nchainlength too high - reducing to %d \n", ULLONG_MAX);
            chainlength = ULONG_MAX;
        }
        else
        {
            chainlength = static_cast<PosLargeInt> (cl);
        }
    }
    else
    {
        // we will do an exhaustive model search, so initialize chainlength
        // with a non-usable number:
        chainlength = 0;
    }
}


// ***************************************************************************************************//

// GlmModelInfo //

// convert to Rcpp list
List
GlmModelInfo::convert2list(long double logNormConst,
                           const Book& bookkeep) const
{
    return List::create(_["logMargLik"] = logMargLik,
                        _["logPrior"] = logPrior,
                        _["posterior"] = NumericVector::create(exp(logPost - logNormConst),
                                                               bookkeep.doSampling ? (hits * 1.0 / bookkeep.chainlength) : NA_REAL),
                        _["negLogUnnormZDensities"] = negLogUnnormZDensities.convert2list(),
                        _["zMode"] = zMode,
                        _["zVar"] = zVar,
                        _["laplaceApprox"] = laplaceApprox,
                        _["residualDeviance"] = residualDeviance);
}


// ***************************************************************************************************//

bool ModelPar::operator<(const ModelPar& m) const
{ 
	// return size() < m.size(); way too easy... lexicographical comparison, starting with uc indices:
	if (ucPars < m.ucPars)
		return true;
	else if (ucPars > m.ucPars)
		return false;
	else // uc indices are equal
		return lexicographical_compare(fpPars.begin(), fpPars.end(), m.fpPars.begin(), m.fpPars.end());	
}

// ***************************************************************************************************//

// return a textual description of this model configuration
std::string
ModelPar::print(const FpInfo& fpInfo) const
{
    // put everything into a stringstream, and return the string conversion at the end.
    std::ostringstream stream;

    // ask for two significant digits
    stream.precision(2);

    // Start
    stream << "\nmodel with " << fpSize << " FP powers and " << ucPars.size() << " UC groups.\n\nIncluded FP powers:";

    if(fpSize > 0)
    {
        // in parallel to the power index vector, go along the names vector
        StrVector::const_iterator n = fpInfo.fpnames.begin();
        for(PowersVector::const_iterator
                i = fpPars.begin();
                i != fpPars.end();
                ++i, ++n)
        {
            // if there is at least one power
            if(i->size() > 0)
            {
                // print the name
                stream << "\nFor " << *n << ": ";

                // and all powers
                MyDoubleVector thesePowers = fpInfo.inds2powers(*i);
                for(MyDoubleVector::const_iterator
                        p = thesePowers.begin();
                        p != thesePowers.end();
                        ++p)
                {
                    stream << *p << " ";
                }
            }
        }
    }
    else
    {
        stream << "\nnone";
    }


    stream << "\n\nIndexes of included UC groups:";

    if(ucPars.size() > 0)
    {
        for(IntSet::const_iterator
                i = ucPars.begin();
                i != ucPars.end();
                ++i)
        {
            stream << " " << *i;
        }
    }
    else
    {
        stream << " none";
    }
    
    stream << "\n\nIndexes of included fixed groups:";
    
    if(fixPars.size() > 0)
    {
      for(IntSet::const_iterator
            i = fixPars.begin();
          i != fixPars.end();
          ++i)
      {
        stream << " " << *i;
      }
    }
    else
    {
      stream << " none";
    }

    // return the resulting string
    return stream.str();
}

// ***************************************************************************************************//

// convert R list to ModelPar
ModelPar::ModelPar(List rcpp_configuration,
                   const FpInfo& fpInfo) :
                   fpSize(0),
                   ucPars(as<IntSet>(rcpp_configuration["ucTerms"])),
                   fixPars(as<IntSet>(rcpp_configuration["fixTerms"]))
{
    // get the list with the FP power vectors:
    List rcpp_powers = rcpp_configuration["powers"];

    // now process all fp names saved in the fpInfo, to construct the correct PowersVector
    for(StrVector::const_iterator s = fpInfo.fpnames.begin();
            s != fpInfo.fpnames.end();
            ++s)
    {
        // get the list element with this name,
        // a double vector,
        // and get the Powers object from that (that is, convert double vector to indexes)
        Powers powers(fpInfo.vec2inds(as<MyDoubleVector>(rcpp_powers[*s])));

        // we have started with fpSize = 0, so we can just add the number of powers now
        fpSize += powers.size();

        // save the powers at the right place
        fpPars.push_back(powers);
    }

    // now we are sure that the PowersVector fpPars is parallel to the fpInfo stuff.
}


// ***************************************************************************************************//

List
ModelPar::convert2list(const FpInfo& currFp) const
{
    // powers
    List powers(fpPars.size());
    powers.names() = currFp.fpnames;

    for (PosInt i = 0; i != fpPars.size(); ++i)
    {
        powers[i] = currFp.inds2powers(fpPars[i]);
    }

    // return with uc settings
    return List::create(_["ucTerms"] = ucPars,
                        _["powers"] = powers,
                        _["fixTerms"] = fixPars);
}


// ***************************************************************************************************//

// compute set of free uc group indices in a model configuration
IntSet
ModelPar::getFreeUcs(const PosIntVector& ucSizes,
                     const PosInt& currDim,
                     const PosInt& maxDim) const
{
    IntSet ret;

    for (PosIntVector::size_type i = 1; i <= ucSizes.size(); i++)
    { // for every uc index
        if ((find(ucPars.begin(), ucPars.end(), i) == ucPars.end())
                && (ucSizes.at(i - 1) <= maxDim - currDim))
        {
            ret.insert(i); // insert if not already in model and enough space in design matrix
        }
    }

    return ret;
}

// ***************************************************************************************************//

// compute set of free cov indices in a model configuration
PosIntSet
ModelPar::getFreeCovs(const FpInfo& currFp,
                      const IntSet& freeUcs,
                      const PosInt& currDim,
                      const PosInt& maxDim) const
{
    PosIntSet ret;

    if (currDim == maxDim)
    {
        return ret;
    }

    for (PosInt i = 0; i != fpPars.size(); i++)
    {
        if (fpPars.at(i).size() < currFp.fpmaxs.at(i))
        {
            ret.insert(i + 1);
        }
    }

    if (! freeUcs.empty())
    {
        ret.insert(fpPars.size() + 1);
    }

    return ret;
}

// ***************************************************************************************************//

// determine set of present cov indices
PosIntSet
ModelPar::getPresentCovs() const
{
    PosIntSet ret;

    for (PosInt i = 0; i != fpPars.size(); i++)
    {
        if (! fpPars.at(i).empty())
        {
            ret.insert(i + 1);
        }
    }

    if (! ucPars.empty())
    {
        ret.insert(fpPars.size() + 1);
    }

    return ret;
}


// ***************************************************************************************************//

// push back index into covGroupWisePosteriors-Array
void
ModelPar::pushInclusionProbs(const FpInfo& fpInfo,
                             const UcInfo& ucInfo,
                             Book& bookkeep) const
{
    for (PosInt i = 0; i != fpInfo.nFps; i++)
    {
        if (! fpPars.at(i).empty())
            bookkeep.covGroupWisePosteriors[i].add(bookkeep.modelCounter);
    }

    for (PosInt i = 1; i <= ucInfo.nUcGroups; i++)
    {
        // search for uc group i
        IntSet::const_iterator ipos = find(ucPars.begin(),
                                           ucPars.end(), i);
        if (ipos != ucPars.end())
        { // if mod.ucPars contains i
            bookkeep.covGroupWisePosteriors[i - 1 + fpInfo.nFps].add(
                                                                     bookkeep.modelCounter);
        }
    }
}



// ***************************************************************************************************//

// we do not need this function because two posteriors will never be exactly equal.
// And if they happen to be equal, we do not care about the number of parameters.

//bool Model::operator<(const Model& m) const  // less
//{
//	double thisLogPost = info.logMargLik + info.logPrior;
//	double mLogPost = m.info.logMargLik + m.info.logPrior;
//	if (thisLogPost < mLogPost)
//		return true;
//	else if (thisLogPost > mLogPost)
//		return false;
//	else  // posteriors are equal, then the parameter makes the decision
//		return m.par < par;
//}

// ***************************************************************************************************//

// return the size of the model (excluding the intercept),
// i.e. the number of coefficients
PosInt 
ModelPar::size(const UcInfo& ucInfo, const FixInfo& fixInfo) const 
{
    // number of FP coefficients is easy
    PosInt ret = fpSize;

    // but be careful for the UC coefficients!
    for(IntSet::const_iterator
            g = ucPars.begin();
            g != ucPars.end();
            ++g)
    {
        ret += ucInfo.ucSizes.at(*g - 1);
    }

    // and for fixed coefficients!
    for(IntSet::const_iterator
          g = fixPars.begin();
        g != fixPars.end();
        ++g)
    {
      ret += fixInfo.fixSizes.at(*g - 1);
    }
    
    // return the total
    return ret;
}

// ***************************************************************************************************//

// convert glm model into list for export to R
List
Model::convert2list(const FpInfo& currFp,
                    long double logNormConst,
                    const Book& bookkeep) const
{
    return List::create(_["configuration"] = par.convert2list(currFp),
                        _["information"] = info.convert2list(logNormConst, bookkeep));
}

// ***************************************************************************************************//

// compute nodes and log weights for given mode and var of target unnormalized density
void
GaussHermite::getNodesAndLogWeights(double mode, double var,
                                    MyDoubleVector& nodes, MyDoubleVector& logWeights) const // output
{
    // logarithm of square root of (2 * var).
    double logSqrt2Var = 0.5 * (M_LN2 + log(var));

    MyDoubleVector::const_iterator t = tVec.begin();
    for(MyDoubleVector::const_iterator
            w = wVec.begin();
            w != wVec.end();
            ++w, ++t)
    {
        nodes.push_back(mode + exp(logSqrt2Var) * (*t));
        logWeights.push_back(log(*w) + (*t) * (*t) + logSqrt2Var);
    }
}

// ***************************************************************************************************//

// GlmModelConfig //

// constructor
GlmModelConfig::GlmModelConfig(List& rcpp_family,
                               double nullModelLogMargLik,
                               double nullModelDeviance,
                               double fixedg,
                               S4& rcpp_gPrior,
                               const AVector& responses,
                               bool debug,
                               bool useFixedc,
                               double empiricalMean,
                               bool empiricalgPrior) :
    dispersions(as<NumericVector>(rcpp_family["dispersions"])),
    weights(as<NumericVector>(rcpp_family["weights"])),
    linPredStart(as<NumericVector>(rcpp_family["linPredStart"])),
    offsets(as<NumericVector>(rcpp_family["offsets"])),
    nullModelLogMargLik(nullModelLogMargLik),
    nullModelDeviance(nullModelDeviance),
    fixedg(fixedg),
    familyString(as<std::string>(rcpp_family["family"])),
    linkString(as<std::string>(rcpp_family["link"])),
    canonicalLink((familyString == "binomial" && linkString == "logit") ||
                      (familyString == "poisson" && linkString == "log")),
    empiricalgPrior(empiricalgPrior)
{
    // and the phi from the R family object
    const double phi = rcpp_family["phi"];

    if (familyString == "binomial")
    {
        distribution = new Binomial(responses,
                                    weights);
    }
    else if (familyString == "gaussian")
    {
        distribution = new Gaussian(responses,
                                    weights,
                                    phi);
    }
    else if (familyString == "poisson")
    {
        distribution = new Poisson(responses,
                                   weights);
    }
    else
    {
        Rf_error("Distribution not implemented");
    }

    if (linkString == "logit")
    {
        link = new LogitLink();
    }
    else if (linkString == "probit")
    {
        link = new ProbitLink();
    }
    else if (linkString == "cloglog")
    {
        link = new CloglogLink();
    }
    else if (linkString == "inverse")
    {
        link = new InverseLink();
    }
    else if (linkString == "log")
    {
        link = new LogLink();
    }
    else if (linkString == "identity")
    {
        link = new IdentityLink();
    }
    else
    {
        Rf_error("Link not implemented!");
    }


    if (useFixedc) {
        
         // from the link and the distribution we can derive the constant factor
         // c = v(h(0)) / h'(0)^2
        double deriv = link->mu_eta(0);
        cfactor = distribution->variance(link->linkinv(0)) / (deriv * deriv);
    }
    else if (!useFixedc) {
      	//Rprintf("used y mean of %f\n", empiricalMean);
        double deriv = link->mu_eta(link->linkfun(empiricalMean));
        cfactor = distribution->variance(link->linkinv(link->linkfun(empiricalMean))) / (deriv * deriv);
        
    }

    if(debug)
    {
       Rprintf("Factor c is %f\n", cfactor);
    }

    // ensure that this is positive
    if(! (cfactor > 0.0))
    {
        Rf_error("cfactor equal to %f, so not positive", cfactor);
    }
    
    //we don't need c if we use the empirical prior
    if(empiricalgPrior) cfactor = 1;

    // finally the g-prior stuff:

    // first get the class name of the S4 g-prior object
    std::string gPriorString = rcpp_gPrior.attr("class");

    // and then depending on the name, initialize our gPrior in C++.
    if (gPriorString == "HypergPrior")
    {
        gPrior = new HypergPrior(as<double>(rcpp_gPrior.slot("a")));
    }
    else if (gPriorString == "InvGammaGPrior")
    {
        gPrior = new InvGammaGPrior(as<double>(rcpp_gPrior.slot("a")),
                                    as<double>(rcpp_gPrior.slot("b")));
    }
    else if (gPriorString == "IncInvGammaGPrior")
    {
        gPrior = new IncInvGammaGPrior(as<double>(rcpp_gPrior.slot("a")),
                                       as<double>(rcpp_gPrior.slot("b")));
    }
    else if (gPriorString == "CustomGPrior")
    {
        gPrior = new CustomGPrior(as<SEXP>(rcpp_gPrior.slot("logDens")));
    }
    else
    {
        Rf_error("g-prior not implemented!");
    }
}


// ***************************************************************************************************//



// dataValues //

// 03/12/2012: add censoring indicator vector

DataValues::DataValues(const AMatrix &x,
                       const AMatrix &xcentered,
                       const AVector &y,
                       const IntVector &censInd,
                       const double &totalNum,
                       const IntSet& fixedCols) :
            design(x),
            centeredDesign(xcentered),
            response(y),
            censInd(censInd),
            nObs(design.n_rows),
            onesVector(arma::ones<AVector>(nObs)),
            totalNumber(static_cast<PosLargeInt> (totalNum)),
            fixedCols(fixedCols)
{
    // and the SST
    AVector centeredResponse = response - arma::mean(response);
    sumOfSquaresTotal = arma::dot(centeredResponse, centeredResponse);
}	



// ***************************************************************************************************//


// ModelCache //

// insert model parameter and corresponding info into cache,
// with caring about the maximum number of elements in the map.
bool
ModelCache::insert(const ModelPar& par, const GlmModelInfo& info)
{
    // first check size of cache
    if(isFull())
    {
        // if we are full, then check if this log posterior is better than
        // the worst cached model, which is pointed to by
        MapType::iterator worstModelIter = *(modelIterSet.begin());

        // the comparison
        if((worstModelIter->second) < info)
        {
            // new model is better than worst model cached.
            // so we delete the worst model from the cache.

            // first from the map
            modelMap.erase(worstModelIter);
            // and then from the set
            modelIterSet.erase(modelIterSet.begin());
        }
        else
        {
            // the new model is not better than the worst model cached,
            // so we do not cache it.
            return false;
        }
    }

    // so now we know that we want to insert the model into the cache,
    // either because the cache was not full or because the new model was better
    // than the worst model cached.

    // -> try inserting into the map:
    std::pair<MapType::iterator, bool> ret = modelMap.insert(MapType::value_type(par, info));

    // if we were successful:
    if(ret.second)
    {
        // then also insert the iterator pointing to the map element into the set.
        modelIterSet.insert(ret.first);

        // return success
        return true;
    }
    else
    {
        return false;
        Rf_error("Should not happen: model already contained in model cache!");
    }
}

// search for the log marginal likelihood of a model config in the map,
// and return NA if not found
GlmModelInfo
ModelCache::getModelInfo(const ModelPar& par) const
{
    // search for the config in the map
    MapType::const_iterator ret = modelMap.find(par);

    // if found, return the log marg lik
    if(ret != modelMap.end())
        return ret->second;
    else
        return GlmModelInfo(R_NaReal, R_NaReal, Cache(), 0.0, 0.0, 0.0, R_NaReal);
}

// increment the sampling frequency for a model configuration
// (of course, if this config is not cached nothing is done)
void
ModelCache::incrementFrequency(const ModelPar& par)
{
    // search for the config in the map
    MapType::iterator ret = modelMap.find(par);

    // if found, increment the hits
    if(ret != modelMap.end())
        ret->second.hits++;
}

// compute the log normalising constant from all cached models
long double
ModelCache::getLogNormConstant() const
{
    // use safe summation
    SafeSum vec;

    // traverse the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // and add all unnormalized log posteriors
        vec.add(m->second.logPost);
    }

    // return the log of the sum of the exp'ed saved elements
    return vec.logSumExp();
}

// compute the inclusion probabilities from all cached models,
// taking the log normalising constant, the number of FPs and the number of UC groups
MyDoubleVector
ModelCache::getInclusionProbs(long double logNormConstant, PosInt nFps, PosInt nUcs) const
{
    // abbreviation
    typedef std::vector<SafeSum> SafeSumVector;
    // allocate vector of safeSum objects for all FPs
    SafeSumVector fps(nFps);

    // and all UC groups
    SafeSumVector ucs(nUcs);

    // now process each model in the cache
    for(MapType::const_iterator
            m = modelMap.begin();
            m != modelMap.end();
            ++m)
    {
        // abbrevs
        const ModelPar& thisPar = m->first;
        const GlmModelInfo& thisInfo = m->second;

        // first process the FPs
        {
        SafeSumVector::iterator s = fps.begin();
        for (PowersVector::const_iterator
                p = thisPar.fpPars.begin();
                p != thisPar.fpPars.end();
                ++p, ++s)
        {
            // is this FP in the model m?
            if (! p->empty())
            {
                // then add the normalized model probability onto his FP stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }

        // then process the UC groups
        {
        SafeSumVector::iterator s = ucs.begin();
        for (PosInt i = 1; i <= nUcs; ++i, ++s)
        {
            // is this UC group in the model m?
            if (thisPar.ucPars.find(i) != thisPar.ucPars.end())
            {
                // then add the normalized model probability onto his UC stack
                s->add(exp(thisInfo.logPost - logNormConstant));
            }
        }
        }
    } // end processing all models in the cache

    // so now we can sum up safesum-wise to the return double vector
    MyDoubleVector ret;

    for(SafeSumVector::iterator
            s = fps.begin();
            s != fps.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    for(SafeSumVector::iterator
            s = ucs.begin();
            s != ucs.end();
            ++s)
    {
        ret.push_back(s->sum());
    }

    return ret;
}

// convert the best nModels from the cache into an R list
List
ModelCache::getListOfBestModels(const FpInfo& fpInfo,
                                long double logNormConst,
                                const Book& bookkeep) const
{
    // allocate the return list
    List ret(std::min(bookkeep.nModels,
                      static_cast<PosInt>(modelIterSet.size())));
    // cast is necessary for gcc-4.2 on Mac on R-forge.

    // process the ordered list of best models from the end (because the set is ordered increasingly)
    PosInt i = 0;
    for(SetType::const_reverse_iterator
            s = modelIterSet.rbegin();
            (i < bookkeep.nModels) && (s != modelIterSet.rend());  // so the return list has min(nModels, modelIterSet.size()) elements.
            ++s, ++i)
    {
        // allocate two-element list in the i-th slot of the return list
        ret[i] = List::create(_["configuration"] = (**s).first.convert2list(fpInfo),
                              _["information"] = (**s).second.convert2list(logNormConst,
                                                                           bookkeep));
    }

    return ret;
}


// End of file.

