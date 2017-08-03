// *      13/07/2015 Replace assert() with Rccp:Stop()


#include <iwls.h>
#include <types.h>
#include <rcppExport.h>
//#include <cassert>
#include <design.h>
#include <stdexcept>
#include <sstream>
#include <linalgInterface.h>

// criterion for comparison of two Column vectors of the same size
// max_j (abs(a_j - b_j) / abs(b_j) + 0.01)
// this is similar to the criterion used by R's glm routine on the deviance scale.
// However, we want to avoid the posterior scale because it would slow down the algorithm...
// (computing this criterion is easier than computing the likelihood * prior)
double
criterion(const AVector& a, const AVector& b)
{
    // check lengths
    //assert(a.n_elem == b.n_elem);
    if(a.n_elem != b.n_elem) Rcpp::stop("iwls.cpp:criterion: a.n_elem != b.n_elem");
    
    // this will be the returned value
    double ret = 0.0;

    // now iterate over the elements
#pragma omp parallel for
    for (PosInt j = 0; j < a.n_elem; ++j)
    {
        double tmp = fabs(a(j) - b(j)) / (fabs(b(j)) + 0.01);

        // note that the "critical" directive is vital here!!
        // Otherwise two executions of the same code can lead to different answers.
#pragma omp critical
        ret = (tmp > ret) ? tmp : ret; /* fmax(ret, tmp); */
    }

    // return the criterion value.
    return ret;
}


// constructor: constructs the Iwls object for given model and data.
Iwls::Iwls(const ModelPar &mod,
           const DataValues& data,
           const FpInfo& fpInfo,
           const UcInfo& ucInfo,
           const FixInfo& fixInfo,
           const GlmModelConfig& config,
           const AVector& linPredStart,
           bool conditional,
           double epsilon,
           bool debug,
           bool tbf) :
           design(getDesignMatrix(mod, data, fpInfo, ucInfo, fixInfo)),
           nCoefs(design.n_cols),
           isNullModel(nCoefs == 1),
           useFixedZ(conditional),
           nObs(design.n_rows),
           // not possible because this could be the null model: designWithoutIntercept(nObs, nCoefs - 1),
           response(data.response),
           config(config),
           invSqrtDispersions(1.0 / arma::sqrt(config.dispersions)),
           unscaledPriorPrec(nCoefs, nCoefs),
           results(linPredStart, nCoefs),
           epsilon(epsilon),
           // verbose(debug),
           tbf(tbf)
{
    // only do additional computations if not the TBF methodology is used
    if(! tbf)
    {
        if(! isNullModel)
        {
          AMatrix scaledDesignWithoutInterceptCrossprod;
          
           if(!config.empiricalgPrior){
            // Scaled design matrix without the intercept.
            // This will be diag(dispersions)^(-1/2) * design[, -1]
            // (attention with 0-based indexing of Armadillo objects!)
            AMatrix scaledDesignWithoutIntercept = arma::diagmat(invSqrtDispersions) * design.cols(1, nCoefs - 1);

            // then the log of the determinant of B'(dispersions)^(-1)B, which is part of the submatrix of R^-1
            // (we know that B'(dispersions)^(-1)B is positive definite, so we do not need to check the sign of the determinant)
            scaledDesignWithoutInterceptCrossprod = arma::trans(scaledDesignWithoutIntercept) * scaledDesignWithoutIntercept;
           }
           
            AMatrix infoMatrix;
            if(config.empiricalgPrior)
            {
              IwlsResults resultsFisher = results;
              
              unscaledPriorPrec.zeros();
              infoMatrix = getInformation(100,
                                          nObs,
                                          invSqrtDispersions,
                                          resultsFisher,
                                          config,
                                          response,
                                          design,
                                          epsilon,
                                          unscaledPriorPrec
              ); //an alternative to using the X'X covariance matrix
              
              scaledDesignWithoutInterceptCrossprod = infoMatrix.submat(1, 1, nCoefs - 1, nCoefs - 1);
              
            }
            
            // input that (the main ingredient) into the unscaled R^-1
            // but first be sure there are zeroes anywhere else:
            unscaledPriorPrec.zeros();
            
            unscaledPriorPrec.submat(1, 1, nCoefs - 1, nCoefs - 1) = scaledDesignWithoutInterceptCrossprod / config.cfactor;

           //  Rcpp::Rcout << "Prior Prec:\n"<< unscaledPriorPrec << std::endl;
            
           // if(config.empiricalgPrior) unscaledPriorPrec.submat(1, 1, nCoefs - 1, nCoefs - 1) = infoMatrix.submat(1, 1, nCoefs - 1, nCoefs - 1);
            

            // now directly use the cholesky routine to avoid copying too much unnecessarily
            int info = potrf(false,
                             scaledDesignWithoutInterceptCrossprod);

            // check that all went well
            if(info != 0)
            {
                std::ostringstream stream;
                stream << "dpotrf(scaledDesignWithoutInterceptCrossprod) got error code " << info << "in Iwls constructor";
                throw std::domain_error(stream.str().c_str());
            }
            // now scaledDesignWithoutInterceptCrossprod contains the cholesky factor!

            // also do not copy the cholesky factor saved now in nonInterceptDesignCrossprod into an extra matrix,
            // but use the Triangular view
            logScaledDesignWithoutInterceptCrossprodDeterminant =
                    2.0 * arma::as_scalar(arma::sum(arma::log(arma::diagvec(scaledDesignWithoutInterceptCrossprod))));
        } else {
            // this is the null model

            // be sure that this is correct:
            unscaledPriorPrec(0, 0) = 0.0;
        }
    }
}

// 03/07/2013: add offsets

// do the Iwls algorithm for a given covariance factor g and start linear predictor linPred,
// until convergence or until the maximum number of iterations is reached.
// so also only one iwls step can be performed with this function.
// Note that the linear predictor is the sum of X^T * beta and the vector of offsets.
// returns the number of iterations.
PosInt
Iwls::startWithLastLinPred(PosInt maxIter,
                              double g)
{
    // initialize iteration counter and stopping criterion
    PosInt iter = 0;
    bool converged = false;

    // do IWLS for at most 30 iterations and unfulfilled stopping criterion
    while ((iter++ < maxIter) && (! converged))
    {
        // compute the pseudo-observations and corresponding sqrt(weights) from the linear predictor
        AVector pseudoObs(nObs);
        AVector sqrtWeights(invSqrtDispersions);

#pragma omp parallel for
        for(PosInt i = 0; i < nObs; ++i)
        {
            double mu = config.link->linkinv(results.linPred(i));
            double dmudEta = config.link->mu_eta(results.linPred(i));

            pseudoObs(i) = results.linPred(i) - config.offsets(i) + (response(i) - mu) / dmudEta;
            sqrtWeights(i) *= dmudEta / sqrt(config.distribution->variance(mu));
        }

        // calculate X'sqrt(W), which is needed twice
        AMatrix XtsqrtW = arma::trans(design) * arma::diagmat(sqrtWeights);

        // calculate the precision matrix Q by doing a rank update:
        // if full Bayes is used, then:
        // Q = tcrossprod(X'sqrt(W)) + 1/g * unscaledPriorPrec
        // if TBF are used, then:
        // Q = tcrossprod(X'sqrt(W))
        double scaleFactor = tbf ? 0.0 : 1.0 / g;
        results.qFactor = unscaledPriorPrec;
        syrk(false,
             false,
             XtsqrtW,
             scaleFactor,
             results.qFactor);

        // decompose into Cholesky factor, Q = LL':
        int info = potrf(false,
                         results.qFactor);

        // check that no error occured
        if(info != 0)
        {
            std::ostringstream stream;
            stream << "Cholesky factorization Q = LL' got error code " << info <<
                    " in IWLS iteration " << iter;
            if(! tbf)
            {
                stream << " for z=" << ::log(g);
            }
            throw std::domain_error(stream.str().c_str());
        }

        // save the old coefficients vector
        AVector coefs_old = results.coefs;

        // the rhs of the equation Q * m = rhs   or    R'R * m = rhs
        pseudoObs = arma::diagmat(sqrtWeights) * pseudoObs;
        results.coefs = XtsqrtW * pseudoObs;
        // note that we have some steps to go until the computation
        // of results.coefs is finished!

        // forward-backward solve LL' * v = rhs
        info = potrs(false,
                     results.qFactor,
                     results.coefs);

        // check that no error occured
        if(info != 0)
        {
            std::ostringstream stream;
            stream << "Forward-backward solve got error code " << info <<
                    " in IWLS iteration " << iter;
            if(! tbf)
            {
                stream << " for z=" << ::log(g);
            }
            throw std::domain_error(stream.str().c_str());
        }

        // the new linear predictor is
        results.linPred = design * results.coefs + config.offsets;

        // compare on the coefficients scale, but not in the first iteration where
        // it is not clear from where coefs_old came. Be safe and always
        // decide for non-convergence in this case.
        converged = (iter > 1) ? (criterion(coefs_old, results.coefs) < epsilon) : false;
    }

    // do not (!)
    // warn if IWLS did not converge within the maximum number of iterations
    // because the maximum number of iterations can be set by user of this function.

    // compute log precision determinant
    results.logPrecisionDeterminant = 2.0 * arma::as_scalar(arma::sum(arma::log(arma::diagvec(results.qFactor))));

    // last but not least return the number of iterations
    return iter;
}

// do the Iwls algorithm for a given covariance factor g and new start linear predictor
// linPredStart.
PosInt
Iwls::startWithNewLinPred(PosInt maxIter,
                          double g,
                          const AVector& linPredStart)
{
    // copy start value into linear predictor of the iwls object
    results.linPred = linPredStart;

    // then run the iwls algorithm
    return startWithLastLinPred(maxIter, g);
}

// 03/07/2013: add offsets

// do the Iwls algorithm for a given covariance factor g and new start coefficients vector
// coefsStart.
PosInt
Iwls::startWithNewCoefs(PosInt maxIter,
                  double g,
                  const AVector& coefsStart)
{
    // start with new linpred deriving from the coefs
    return startWithNewLinPred(maxIter, g, design * coefsStart + config.offsets);
}

// 03/07/2013: add offsets

// compute the log of the (unnormalized)
// posterior density for a given parameter consisting of the coefficients vector and z.
//
// Note that it is important to incorporate all model-depending constants here,
// because this function is also used in the Chib-Jeliazkov marginal likelihood estimation,
// comparing different models!!
//
// useFixedZ: is the log-covariance factor z fixed? Then the conditional posterior
// density of the coefficients vector is returned (so the prior of z is not included).
double
Iwls::computeLogUnPosteriorDens(const Parameter& sample) const
{
    // compute the sample of the linear predictor:
    AVector linPredSample = design * sample.coefs + config.offsets;

    // compute the resulting mean vector from the linear predictor via the response function
    AVector meansSample(linPredSample.n_elem);

#pragma omp parallel for
    for(PosInt i = 0; i < meansSample.n_elem; ++i)
    {
        meansSample(i) = config.link->linkinv(linPredSample(i));
    }

    // start with the log likelihood of this coefficients, it is always included
    // this part is included in both cases because it does not depend on
    // the prior on the (non-intercept) effects:

    double ret = config.distribution->loglik(meansSample.memptr());

    // now it depends again on null model or not.

    if(! isNullModel)
    {
        // map z sample on the original g scale
        double g = exp(sample.z);

        // calculate ||(dispersions)^(-1/2) * B * beta||^2
        AVector scaledBcoefsSample = arma::diagmat(invSqrtDispersions) * (linPredSample - sample.coefs(0));

        // this avoids this multiplication of general matrices:
        // "DEVector scaledBcoefsSample = scaledDesignWithoutIntercept * sample.coefs(_(2, nCoefs));"
        double scaledBcoefsSampleNormSquared = arma::dot(scaledBcoefsSample, scaledBcoefsSample);

        // now add the non-null model specific part, which comes from the prior on
        // the coefficients
        ret += 0.5 * (logScaledDesignWithoutInterceptCrossprodDeterminant -
                      scaledBcoefsSampleNormSquared / (g * config.cfactor) -
                     (nCoefs - 1.0) * (2.0 * M_LN_SQRT_2PI + sample.z + log(config.cfactor)));

        if(! useFixedZ)
        {
            // and the log prior of this g
            double logGPrior = config.gPrior->logDens(g);

            // add the log prior of z
            ret += logGPrior + sample.z;

            // note that the sample.z has its origin in the density transformation from g to z.
            // if the prior on g was discrete and not continuous, this part would not be included in general.
        }
    }

    // return the correct value
    return ret;
}

// compute the deviance of the current model, which in R is done by the glm.fit function
// This is required to compute the TBF.
// Note that this changes the Iwls object, because it iterates until convergence
// or until the maximum number of iterations has been reached.
double
Iwls::computeDeviance(PosInt maxIter)
{
    PosInt iter = startWithLastLinPred(maxIter, R_PosInf);

    // if IWLS did not converge within the maximum number of iterations,
    // throw error which is then caught by getGlmVarLogMargLik later on
    if(iter > maxIter)
    {
        std::ostringstream stream;
        stream << "model could not be fitted within " << maxIter <<
                " IWLS iterations. This could be a result of separability in logistic regression, e.g.";
        throw std::domain_error(stream.str().c_str());
    }

    // compute the resulting mean vector from the linear predictor via the response function
    AVector means(nObs);

#pragma omp parallel for
    for(PosInt i = 0; i < nObs; ++i)
    {
        means(i) = config.link->linkinv(results.linPred(i));
    }

    // compute the log-likelihood
    double logLik = config.distribution->loglik(means.memptr());

    // return the deviance, which is just the scaled log-likelihood
    double ret = - 2.0 * logLik;
    return ret;
}


// Compute a standard GLM to get the observed Fisher Information to use as covariance matrix.
AMatrix 
Iwls::getInformation(PosInt maxIter, PosInt nObs, AVector invSqrtDispersions, IwlsResults results,
          const GlmModelConfig& config, const AVector& response, const AMatrix design,
          double epsilon, AMatrix unscaledPriorPrec)
{
  
  // initialize iteration counter and stopping criterion
  PosInt iter = 0;
  bool converged = false;
  
  AVector finalWeights(invSqrtDispersions);
  
  // do IWLS for at most 30 iterations and unfulfilled stopping criterion
  while ((iter++ < maxIter) && (! converged))
  {
    // compute the pseudo-observations and corresponding sqrt(weights) from the linear predictor
    AVector pseudoObs(nObs);
    AVector sqrtWeights(invSqrtDispersions);
    
#pragma omp parallel for
    for(PosInt i = 0; i < nObs; ++i)
    {
      double mu = config.link->linkinv(results.linPred(i));
      double dmudEta = config.link->mu_eta(results.linPred(i));
      
      pseudoObs(i) = results.linPred(i) - config.offsets(i) + (response(i) - mu) / dmudEta;
      sqrtWeights(i) *= dmudEta / sqrt(config.distribution->variance(mu));
    }
    
    // calculate X'sqrt(W), which is needed twice
    AMatrix XtsqrtW = arma::trans(design) * arma::diagmat(sqrtWeights);
    
    // calculate the precision matrix Q by doing a rank update:
    // if full Bayes is used, then:
    // Q = tcrossprod(X'sqrt(W)) + 1/g * unscaledPriorPrec
    // if TBF are used, then:
    // Q = tcrossprod(X'sqrt(W))
    double scaleFactor = 0.0;
    results.qFactor = unscaledPriorPrec; //just zeros to set up matrix
    syrk(false,
         false,
         XtsqrtW,
         scaleFactor,
         results.qFactor);
    
    // decompose into Cholesky factor, Q = LL':
    int info = potrf(false,
                     results.qFactor);
    
    // check that no error occured
    if(info != 0)
    {
      std::ostringstream stream;
      stream << "Cholesky factorization Q = LL' got error code " << info <<
        " in IWLS iteration " << iter;

      throw std::domain_error(stream.str().c_str());
    }
    
    // save the old coefficients vector
    AVector coefs_old = results.coefs;
    
    // the rhs of the equation Q * m = rhs   or    R'R * m = rhs
    pseudoObs = arma::diagmat(sqrtWeights) * pseudoObs;
    results.coefs = XtsqrtW * pseudoObs;
    
    // note that we have some steps to go until the computation
    // of results.coefs is finished!
    
    // forward-backward solve LL' * v = rhs
    info = potrs(false,
                 results.qFactor,
                 results.coefs);
    
    // check that no error occured
    if(info != 0)
    {
      std::ostringstream stream;
      stream << "Forward-backward solve got error code " << info <<
        " in IWLS iteration " << iter;

      throw std::domain_error(stream.str().c_str());
    }
    
    // the new linear predictor is
    results.linPred = design * results.coefs + config.offsets;
    
    // compare on the coefficients scale, but not in the first iteration where
    // it is not clear from where coefs_old came. Be safe and always
    // decide for non-convergence in this case.
    converged = (iter > 1) ? (criterion(coefs_old, results.coefs) < epsilon) : false;
    finalWeights = sqrtWeights;
  }
  
  // do not (!)
  // warn if IWLS did not converge within the maximum number of iterations
  // because the maximum number of iterations can be set by user of this function.
  
  // compute log precision determinant
  results.logPrecisionDeterminant = 2.0 * arma::as_scalar(arma::sum(arma::log(arma::diagvec(results.qFactor))));
  
  
  //to calculate the observed Information we need the dispersion
  double dispersion;
  if( config.familyString == "poisson" || config.familyString == "binomial" ){
    dispersion = 1;
  } else {
    dispersion = dot(finalWeights, arma::square(response - results.linPred))/(nObs-results.coefs.n_elem);
  }
  
  
  AMatrix observed = arma::trans(design) *  arma::diagmat(arma::square(finalWeights)) * design / dispersion;
  
  
  AVector residuals = arma::square(response - results.linPred);
  
  return(observed);
}
