#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <set>
#include <map>
#include <vector>
#include <iterator>
#include <numeric>
#include <string>

#include <rcppExport.h>
#include <types.h>
#include <functionWraps.h>
#include <fpUcHandling.h>
#include <links.h>
#include <distributions.h>
#include <gpriors.h>


// ***************************************************************************************************//

struct SafeSum
{
    LongDoubleVector vals;

    void
    add(const long double &val);

    // compute the sum of the elements using accurate algorithm
    long double
    sum();

    // compute the log of the sum of the exp of the elements using accurate algorithm,
    // and avoiding infinite contributions.
    long double
    logSumExp();

    // compute the sum of the elements using very simple algorithm
    long double
    simpleSum();
};

// ***************************************************************************************************//

struct IndexSafeSum
{
    // the model index collection
    IndSet indices;

    // add a model index i to the collection indices
    void
    add(const Ind& i);

    // taking the safe sum object s of log posteriors, and the associated log normalizing
    // constant logNormConst,
    // compute the sum \sum_i in indices exp{s_i - logNormConst}.
    long double
    sumNormalizedExp(const SafeSum& s, long double logNormConst) const;
};

// ***************************************************************************************************//

// 21/11/2012: add tbf option

// bookkeeping, shall later be used by glm and hyper-g!
struct Book
{
    PosLargeInt modelCounter;
    PosLargeInt chainlength;
    PosLargeInt nanCounter;

    SafeSum modelLogPosteriors;
    //IndexSafeSum *covGroupWisePosteriors; // for computation of covariate inclusion probs: array (bfp, uc)
    std::vector<IndexSafeSum> covGroupWisePosteriors; // for computation of covariate inclusion probs: array (bfp, uc)

    const bool tbf;
    const bool doGlm;
    const bool empiricalBayes;
    const bool useFixedg;
    const bool useFixedc;
    const bool doSampling;
    const bool verbose;
    const std::string modelPrior;

    const PosInt nModels;
    const PosInt nCache;

    const double largeVariance;
    const bool useBfgs;
    const bool debug;

    const bool higherOrderCorrection;

    // constructor which checks the chainlength
    Book(bool tbf,
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
         bool higherOrderCorrection);

};

// ***************************************************************************************************//

// data structure: model information for general models.
// keep in mind that contents must be assignable.
struct ModelInfo
{
    // most important:
    double logMargLik;
    double logPrior;

    // just the sum of logMargLik and logPrior.
    // But we need it often, so save it here (the ctr computes it).
    double logPost;

    // only needed for MCMC:
    PosLargeInt hits;

    // simple constructor
    ModelInfo(double logMargLik, double logPrior, PosLargeInt hits = 0) :
        logMargLik(logMargLik), logPrior(logPrior), logPost(logMargLik + logPrior), hits(hits)
    {
    }

    // compare two model infos bei their posterior probability:
    // the model info with the lower one is lower.
    bool
    operator<(const ModelInfo& m) const
    {
        return logPost < m.logPost;
    }
};

// ***************************************************************************************************//

//// specialization for normal model with hyper-g prior
//struct NormalModelInfo : public ModelInfo
//{
//    // additional to ModelInfo:
//
//    // posterior expected factor g given this model
//    double postExpectedg;
//
//    // posterior expected shrinkage factor g/(1+g) given this model
//    double postExpectedShrinkage;
//
//    // coefficient of determination for this model
//    double R2;
//
//    // constructor
//    NormalModelInfo(double logMargLik,
//                    double logPrior,
//                    double postExpectedg,
//                    double postExpectedShrinkage,
//                    double R2,
//                    PosLargeInt hits) :
//        ModelInfo(logMargLik, logPrior, hits),
//        postExpectedg(postExpectedg),
//        postExpectedShrinkage(postExpectedShrinkage),
//        R2(R2)
//        {
//        }
//};

// ***************************************************************************************************//

// specialization for GLM model
struct GlmModelInfo : public ModelInfo
{
    GlmModelInfo(double logMargLik,
                 double logPrior,
                 const Cache& cache,
                 double zMode,
                 double zVar,
                 double laplaceApprox,
                 double residualDeviance) :
                     ModelInfo(logMargLik,
                               logPrior),
                     negLogUnnormZDensities(cache),
                     zMode(zMode),
                     zVar(zVar),
                     laplaceApprox(laplaceApprox),
                     residualDeviance(residualDeviance)
                     {
                     }

    // cache of negative log unnormalized posterior density of z evaluations:

    // the z values
    Cache negLogUnnormZDensities;

    // moment estimates:

    // mode of the z posterior (z*)
    double zMode;

    // and associated estimated variance (sigma*^2)
    double zVar;

    // here is the resulting Laplace approximation of the log marginal likelihood
    double laplaceApprox;

    // this is only filled in the TBF case
    double residualDeviance;

    // convert to R list
    Rcpp::List
    convert2list(long double logNormConst,
                 const Book& bookkeep) const;

    // conversion of R list to modelInfo:
    // Cache and hits are not filled!
    explicit
    GlmModelInfo(Rcpp::List rcpp_information) :
        ModelInfo(rcpp_information["logMargLik"],
                  rcpp_information["logPrior"],
                  NA_INTEGER),
                  zMode(rcpp_information["zMode"]),
                  zVar(rcpp_information["zVar"]),
                  laplaceApprox(rcpp_information["laplaceApprox"]),
                  residualDeviance(rcpp_information["residualDeviance"])
    {
    }
};

// ***************************************************************************************************//

// model parameter object: must have a strict weak ordering
struct ModelPar
{
    PowersVector fpPars; // vector of multisets

    // not needed: just use fpPars.size()
    // PosInt nFps; // length of vector

    PosInt fpSize; // number of fp powers
    IntSet ucPars; // set of group indices, starting from 1 (!)
    IntSet fixPars; // set of group indices, starting from 1 (!)

    // not needed: just use ucPars.size()
    // PosInt ucSize; // number of uc Groups included

    // start with an empty (null) model:
    ModelPar(PosInt nFps) :
        fpPars(nFps), fpSize(0), ucPars(), fixPars()
    {
    }

    // comparison of configurations
    bool
    operator<(const ModelPar& m) const;

    // return a textual description of this model configuration
    std::string
    print(const FpInfo& fpInfo) const;

    // return the size of the model (excluding the intercept)
    PosInt
    size(const UcInfo& ucInfo, const FixInfo& fixInfo) const;

    // convert to R list
    Rcpp::List
    convert2list(const FpInfo& currFp) const;

    // convert R list to ModelPar
    ModelPar(Rcpp::List rcpp_configuration,
             const FpInfo& fpInfo);

    // compute set of free uc group indices in a model configuration
    IntSet
    getFreeUcs(const PosIntVector& ucSizes,
               const PosInt& currDim,
               const PosInt& maxDim) const;

    // compute set of free cov indices in a model configuration
    PosIntSet
    getFreeCovs(const FpInfo& currFp,
                const IntSet& freeUcs,
                const PosInt& currDim,
                const PosInt& maxDim) const;

    // determine set of present cov indices
    PosIntSet
    getPresentCovs() const;

    // push back index into covGroupWisePosteriors-Array
    void
    pushInclusionProbs(const FpInfo& fpInfo,
                       const UcInfo& ucInfo,
                       Book& bookkeep) const;
};

// ***************************************************************************************************//


// first only for Glms. Possibly extended with template class for ModelInfo to NormalModelInfo.
struct Model
{
    ModelPar par;
    GlmModelInfo info;

    // initialize
    Model(const ModelPar& p, const GlmModelInfo& i) :
        par(p), info(i)
    {
    }

    // fully forward the comparison to the model parameter and info class.
    // we need to include the par comparison because sometimes
    // really the info is identical...
    bool
    operator<(const Model& m) const // less
    {
        if(info < m.info)
            return true;
        else if(m.info < info)
            return false;
        else
            return par < m.par;
    }

    // model to list
    Rcpp::List
    convert2list(const FpInfo& currFp,
                 long double logNormConst,
                 const Book& bookkeep) const;
};


// ***************************************************************************************************//


class GaussHermite {
public:
    // compute nodes and log weights for given mode and sd of target unnormalized density
    void
    getNodesAndLogWeights(double mode, double var,
                       MyDoubleVector& nodes, MyDoubleVector& logWeights) const; // output

    // constructor
    GaussHermite(Rcpp::List rcpp_gaussHermite) :
        tVec(Rcpp::as<MyDoubleVector>(rcpp_gaussHermite["nodes"])),
        wVec(Rcpp::as<MyDoubleVector>(rcpp_gaussHermite["weights"]))
        {
        }

private:
    const MyDoubleVector tVec; // nodes
    const MyDoubleVector wVec; // weights (not log!)
};


// ***************************************************************************************************//



struct GlmModelConfig
{
    // vector of dispersions (phi / weights)
    const AVector dispersions;

    // vector of weights
    const AVector weights;

    // vector of starting values for the linear predictor
    // constant, because for each model we will initially start from that.
    const AVector linPredStart;

    // the vector of offsets
    const AVector offsets;

    // log marg lik of the the null model
    const double nullModelLogMargLik;

    // deviance of the the null model
    const double nullModelDeviance;

    // fixed value of g
    const double fixedg;

    // the g-prior information
    const GPrior* gPrior;

    // the link information
    const Link* link;

    // the distribution information
    const Distribution* distribution;

    // the constant factor deriving from the link and the distribution,
    // which we need for the generalized g-prior
    double cfactor;

    // save also the names of link and distribution
    const std::string familyString;
    const std::string linkString;

    // does this model use the canonical link?
    const bool canonicalLink;

    // should we use the empirical g prior which uses 
    // the information matrix ( J(B)^-1 ) for the covariance instead of (X'X)^-1
    const bool empiricalgPrior;
    
    // constructor
    GlmModelConfig(Rcpp::List& rcpp_family,
                   double nullModelLogMargLik,
                   double nullModelDeviance,
                   double fixedg,
                   Rcpp::S4& rcpp_gPrior,
                   const AVector& responses,
                   bool debug,
                   bool useFixedc,
                   double empiricalMean,
                   bool empiricalgPrior);

    // destructor
    ~GlmModelConfig()
    {
        delete gPrior;
        delete link;
        delete distribution;
    }
};

// ***************************************************************************************************//


struct DataValues
{
    AMatrix design;
    AMatrix centeredDesign;

    AVector response;
    double sumOfSquaresTotal;

    IntVector censInd;

    int nObs;

    AVector onesVector;

    PosLargeInt totalNumber; // cardinality of model space

    const IntSet fixedCols;

    DataValues(const AMatrix &x,
               const AMatrix &xcentered,
               const AVector &y,
               const IntVector &censInd,
               const double &totalNum,
               const IntSet& fixedCols);
};


// ***************************************************************************************************//

// first only for GLM models:
// the model cache class.

// Caches the best models in a map of a given maximum size, and also stores the
// (unnormalized) log posterior probabilities in an ordered set, pointing to the models in the map.
class ModelCache {
public:


    // create a new ModelCache with given maximum size.
    ModelCache(int maxSize) :
        maxSize(maxSize),
        modelMap(),
        modelIterSet()
        {
        }

    // check if max size was reached
    bool
    isFull() const
    {
        return modelMap.size() == maxSize;
    }

    // return size of cache
    int
    size() const
    {
        return modelMap.size();
    }

    // insert model parameter and belonging model info into the cache.
    // returns false if not inserted (e.g. because the par was
    // already inside, or the model was not good enough)
    bool
    insert(const ModelPar& par, const GlmModelInfo& info);

    // search for the model info of a model config in the map,
    // and return an information with NA for log marg lik if not found
    GlmModelInfo
    getModelInfo(const ModelPar& par) const;

    // increment the sampling frequency for a model configuration
    // (of course, if this config is not cached nothing is done!)
    void
    incrementFrequency(const ModelPar& par);

    // compute the log normalising constant from all cached models
    long double
    getLogNormConstant() const;

    // compute the inclusion probabilities from all cached models,
    // taking the log normalising constant and the total number of FPs / UC groups
    MyDoubleVector
    getInclusionProbs(long double logNormConstant, PosInt nFps, PosInt nUcs) const;

    // convert the best nModels from the cache into an R list
    Rcpp::List
    getListOfBestModels(const FpInfo& fpInfo,
                        long double logNormConst,
                        const Book& bookkeep) const;


private:
    // the map type
    typedef std::map<ModelPar, GlmModelInfo> MapType;

    // define comparison function for iterators
    struct Compare_map_iterators
    {
        bool
        operator()(const MapType::iterator& first, const MapType::iterator& second) const
        {
            return (first->second.logPost) < (second->second.logPost);
        }
    };

    // the set type of ordered map iterators
    typedef std::set<MapType::iterator, Compare_map_iterators> SetType;

    // and finally the data members
    const MapType::size_type maxSize;
    MapType modelMap;
    SetType modelIterSet;
};

// ***************************************************************************************************//

// all information needed in mcmc function
struct ModelMcmc
{
    // initialize with the null model only.
    ModelMcmc(const FpInfo& fpInfo,
              const UcInfo& ucInfo,
              PosInt maxDim,
              double logMargLikNullModel) :
                  modPar(fpInfo.nFps),
                  dim(1),
                  freeUcs(modPar.getFreeUcs(ucInfo.ucSizes, dim, maxDim)),
                  freeCovs(modPar.getFreeCovs(fpInfo, freeUcs, dim, maxDim)),
                  presentCovs(modPar.getPresentCovs()),
                  birthprob(1.0),
                  deathprob(0.0),
                  moveprob(0.0),
                  logMargLik(logMargLikNullModel),
                  logPrior(R_NaN)
    {
    }

    ModelPar modPar;
    PosInt dim; // number of columns in this model's design matrix

    IntSet freeUcs; // indices within uc groups, denoting the birthable ones
    PosIntSet freeCovs; // indices of free covs (starting from first fp with index 1 up to uc index = nFps + 1)
    PosIntSet presentCovs; // analogue

    double birthprob, deathprob, moveprob; // move type probabilites, switchprob is 1-bprob-dprob-mprob.
    double logMargLik;
    double logPrior;
};

// ***************************************************************************************************//

#endif /*DATASTRUCTURE_H_*/
