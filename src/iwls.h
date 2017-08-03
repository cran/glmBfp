/*
 * zdensity.h
 *
 *  Created on: 16.11.2009
 *      Author: daniel
 */

#ifndef IWLS_H_
#define IWLS_H_

#include <design.h>
#include <types.h>


// only one parameter set for the GLM
struct Parameter
{
    // ctr
    Parameter(const AVector& coefs,
              double z) :
                  coefs(coefs),
                  z(z)
                  {
                  }

    // default ctr
    Parameter(PosInt nCoefs) :
        coefs(nCoefs),
        z(R_NaReal)
        {
        }

    // the sampled regression coefficients (including the intercept)
    AVector coefs;

    // the sampled log covariance factor
    double z;
};




// all info from the iwls results
struct IwlsResults
{
    // ctr
    IwlsResults(const AVector& linPred,
                PosInt nCoefs) :
        linPred(linPred),
        coefs(nCoefs),
        qFactor(nCoefs, nCoefs)
        {
        }

    // default ctr
    IwlsResults(PosInt nObs,
                PosInt nCoefs) :
        linPred(nObs),
        coefs(nCoefs),
        qFactor(nCoefs, nCoefs)
        {
        }

    // the linear predictor (including offsets!)
    AVector linPred;

    // the corresponding coefficient vector,
    AVector coefs;

    // lower-triangular Cholesky factor of the precision matrix
    AMatrix qFactor;

    // the log determinant of the precision matrix
    double logPrecisionDeterminant;
};


// 03/07/2013: add offsets

class Iwls {

public:
    // constructor: constructs the Iwls object for given model and data and
    // start linear predictor
    Iwls(const ModelPar &mod,
         const DataValues& data,
         const FpInfo& fpInfo,
         const UcInfo& ucInfo,
         const FixInfo& fixInfo,
         const GlmModelConfig& config,
         const AVector& linPredStart,
         bool conditional,
         double epsilon,
         bool debug,
         bool tbf);

    // do the Iwls algorithm for a given covariance factor g and current linear predictor
    // linPred,
    // until convergence or until the maximum number of iterations is reached.
    // so also only one iwls step can be performed with this function.
    // Note that the linear predictor is the sum of X^T * beta and the vector of offsets.
    // returns the number of iterations.
    PosInt
    startWithLastLinPred(PosInt maxIter,
                         double g);

    // do the Iwls algorithm for a given covariance factor g and new start linear predictor
    // linPredStart.
    PosInt
    startWithNewLinPred(PosInt maxIter,
                        double g,
                        const AVector& linPredStart);

    // do the Iwls algorithm for a given covariance factor g and new start coefficients vector
    // coefsStart.
    PosInt
    startWithNewCoefs(PosInt maxIter,
                      double g,
                      const AVector& coefsStart);


    // compute the log of the (unnormalized)
    // posterior density for a given parameter consisting of the coefficients vector and z
    double
    computeLogUnPosteriorDens(const Parameter& sample) const;

    // compute the deviance of the current model, which in R is done by the glm.fit function
    // This is required to compute the TBF.
    // Note that this changes the Iwls object, because it iterates until convergence
    // or until the maximum number of iterations has been reached.
    double
    computeDeviance(PosInt maxIter);

    // getter for results
    IwlsResults
    getResults() const
    {
        return results;
    }

    // Get the fisher information for the desired model
    AMatrix getInformation(PosInt maxIter,
                           PosInt nObs,
                           AVector invSqrtDispersions,
                           IwlsResults results,
                           const GlmModelConfig& config,
                           const AVector& response,
                           const AMatrix design,
                           double epsilon,
                           AMatrix unscaledPriorPrec
    );
    
    // this can be public:

    // design matrix for this model
    const AMatrix design;

    // dimension including the intercept ("ncol(design)" in R syntax)
    const PosInt nCoefs;

    // is this the null model?
    const bool isNullModel;

    // use a fixed z?
    const bool useFixedZ;

    // number of observations
    const PosInt nObs;


private:


    // the log of the determinant of the crossproduct of scaledDesignWithoutIntercept
    // This is B'(dispersions)^(-1)B.
    double logScaledDesignWithoutInterceptCrossprodDeterminant;

    // the response vector
    const AVector& response;

    // also needed in function calls (mean, variance and other glm functions)
    const GlmModelConfig& config;

    // This is diag(dispersions)^(-1/2), use diagmat to interpret it as diagonal matrix!
    AVector invSqrtDispersions;

    // unscaled prior precision matrix R^-1 = blockDiag(0, 1/c B'(dispersions)^(-1)B) (without g^-1)
    AMatrix unscaledPriorPrec;

    // container for the results of the iwls computation:
    IwlsResults results;

    // the convergence epsilon
    const double epsilon;

    // status messages??
    // const bool verbose;

    // use TBF methodology?
    const bool tbf;
};


#endif /* IWLS_H_ */
