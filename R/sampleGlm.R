#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[sampleGlm.R] by DSB Mon 26/08/2013 15:43 (CEST)>
##
## Description:
## Produce posterior samples from one GLM returned by glmBayesMfp, using an MCMC sampler. 
##
## History:
## 10/12/2009   file creation
## 05/01/2010   construct the random number generator for the marginal z approximation
##              *inside* the helper function, so that we can use the normalizing constant 
##              computed by Runuran instead of doing another "integrate" call.
## 06/01/2010   add postprocessing of the coefficients samples
## 16/02/2010   correct title of getLogMargLikEstimate man page
## 17/02/2010   split off getMarginalZ into separate code file
## 18/02/2010   split off getLogMargLikEstimate into separate code file
## 15/03/2010   add check for class of mcmc argument ...
## 17/03/2010   ... but do it correctly (the class is called McmcOptions and not Mcmc)
## 08/04/2010   Add description of the return list (S3 class "GlmBayesMfpSamples").
## 17/05/2010   Add "useOpenMP" option, which makes it easier
##              to switch than setting the environment variable OMP_NUM_THREADS
##              (especially when using R on the server via Emacs Tramp...)
## 24/05/2010   return samples on the response scale ("response") instead of their
##              averages ("fitted"), because this averaging can better be done in
##              a "fitted" method and we need the raw samples in sampleBma.R
## 25/05/2010   different structure of return list;
##              layout of bfp and uc samples is now nPar x nSamples.
## 28/05/2010   - remove argument "weights" which is no longer needed here,
##              as we only sample the linear predictors for new observations.
##              - rename "response" to "fitted", which shall also contain the
##              *linear predictors* for the fitted data (and not the *means*).
## 02/07/2010   extend the functionality to the case with fixed z (option "fixedZ")
## 08/07/2010   check if the models were fit by empirical Bayes. In that
##              case use the fixed z option !
## 26/07/2010   Do not drop dimensions when selecting UC coefficients samples
## 23/11/2012   modifications to accommodate the TBF methodology
## 03/12/2012   modifications to accommodate the Cox models
## 24/01/2013   adapt for fixedg option
## 03/07/2013   comment on offsets
#####################################################################################

##' @include helpers.R
##' @include McmcOptions-class.R
##' @include McmcOptions-methods.R
##' @include GlmBayesMfpSamples-class.R
##' @include getMarginalZ.R
##' @include getLogMargLikEstimate.R
{}

##' Produce posterior samples from one GLM / Cox model
##'
##' Based on the result list from \code{\link{glmBayesMfp}}, for the first model
##' in the list MCMC samples are produced. In parallel to the sampling of
##' coefficients and FP curve points, optionally the marginal likelihood of the
##' model is estimated with MCMC samples. This provides a check of the
##' integrated Laplace approximation used in the model sampling. If TBF
##' methodology is used, then no MCMC is necessary, instead ordinary Monte Carlo
##' samples from an approximate posterior distribution are obtained.
##'
##' @param object the \code{GlmBayesMfp} object, from which only the first model
##' will be processed (at least for now \ldots)
##' @param mcmc MCMC options object with class \code{\linkS4class{McmcOptions}}.
##' If TBF is used, each sample is accepted, and the number of samples is given
##' by \code{\link{sampleSize}}(\code{mcmc}).
##' @param estimateMargLik shall the marginal likelihood be estimated in
##' parallel? (default) Only has an effect if full Bayes and not TBF is used. 
##' @param gridList optional list of appropriately named grid vectors for FP
##' evaluation. Default is length (\code{gridSize} - 2) grid per covariate
##' additional to the observed values (two are at the endpoints) 
##' @param gridSize see above (default: 203)
##' @param newdata new covariate data.frame with exactly the names (and
##' preferably ranges) as before (default: no new covariate data) Note that
##' there is no option for offsets for new data at the moment. Just add the
##' offsets to the \code{predictions} slot of \code{samples} in the return list 
##' yourself.
##' @param fixedZ either \code{NULL} (default) or a (single) fixed z value to
##' be used, in order to sample from the conditional posterior given this z.
##' If \code{object} was constructed by the empirical Bayes machinery,
##' this will default to the estimated z with maximum conditional marginal
##' likelihood. If \code{object} was constructed with the option \code{fixedg},
##' then the fixed value will be used by default.
##' @param marginalZApprox method for approximating the marginal density of the
##' log covariance factor z, see \code{\link{getMarginalZ}} for the details
##' (default: same preference list as in \code{\link{getMarginalZ}})
##' If TBF are used in conjunction with incomplete inverse gamma hyperprior on
##' g = exp(z), then the posterior distribution of g is again of this form.
##' Therefore this option does not have any effect in that case, because the
##' samples are directly obtained from that posterior distribution.
##' @param verbose should information on computation progress be given?
##' (default)
##' @param debug print debugging information? (not default)
##' @param useOpenMP shall OpenMP be used to accelerate the computations?
##' (default)
##' @param correctedCenter If TRUE predict new data based on the centering 
##' of the original data.
##' 
##' @return Returns a list with the following elements:
##' \describe{
##' \item{samples}{an object of S4 class
##' \code{\linkS4class{GlmBayesMfpSamples}}} 
##' \item{coefficients}{samples of all original coefficients in the model
##' (nCoefs x nSamples)}
##' \item{acceptanceRatio}{proportion of accepted Metropolis-Hastings proposals}
##' \item{logMargLik}{if \code{estimateMargLik} is \code{TRUE}, this list is
##' included: it contains the elements \code{numeratorTerms} and
##' \code{denominatorTerms} for the numerator and denominator samples of the
##' Chib Jeliazkov marginal likelihood estimate,
##' \code{highDensityPointLogUnPosterior} is the log unnormalized posterior
##' density at the fixed parameter and the resulting \code{estimate} and
##' \code{standardError}.}
##' }
##'
##' @importFrom Rcpp evalCpp
##'
##' @export
##' @keywords models regression
sampleGlm <-
    function(object,
             mcmc=McmcOptions(),
             estimateMargLik=TRUE,
             gridList=list(),
             gridSize=203L,
             newdata=NULL,
             fixedZ=NULL,
             marginalZApprox=NULL,
             verbose=TRUE,
             debug=FALSE,
             useOpenMP=TRUE,
             correctedCenter=FALSE)
{
    ## check the object
    if(! inherits(object, "GlmBayesMfp"))
        stop(simpleError("object must be of class GlmBayesMfp"))
    if(! (length(object) >= 1L))
        stop(simpleError("there has to be at least one model in the object"))

    ## other checks
    stopifnot(is.bool(estimateMargLik),
              is.bool(verbose),
              is.bool(debug),
              is(mcmc, "McmcOptions"),
              is.bool(useOpenMP))
    
    ## coerce newdata to data frame
    newdata <- as.data.frame(newdata)
    nNewObs <- nrow(newdata)

    ## ## covariates matrix for newdata (if there is any):
    ## newX <-
    ##     if(nNewObs > 0L)
    ##     {
    ##         constructNewdataMatrix(object=object,
    ##                                newdata=newdata)
    ##     }
    ##     else
    ##         NULL

    ## get the old attributes of the object
    attrs <- attributes(object)

    ## is this a GLM or a Cox model?
    doGlm <- attrs$distribution$doGlm

    ## take only the first model
    model <- object[[1L]]
    config <- model$configuration
    info <- model$information

    ## get the (centered) design matrix of the model for the original data
    design <- getDesignMatrix(modelConfig=config,
                              object=object,
                              intercept=doGlm) 
    ## the model dimension (including the intercept for GLMs)
    modelDim <- ncol(design)
    
    ## check if the model is the null model
    isNullModel <- identical(modelDim,
                             ifelse(doGlm,
                                    1L,
                                    0L))

    ## check / modify the fixedZ option
    if(isNullModel &&
       is.null(fixedZ))
    {
        fixedZ <- 0
    }
    else if((attrs$searchConfig$useFixedg ||
             attrs$searchConfig$empiricalBayes) &&
            is.null(fixedZ))
    {
        fixedZ <- info$zMode
    }
    
    if(useFixedZ <- ! is.null(fixedZ))
    {
        stopifnot(is.numeric(fixedZ),
                  is.finite(fixedZ),
                  identical(length(fixedZ), 1L))
    } else {
        ## give some value to fixedZ for the C++ code side...
        fixedZ <- 0
    }

    ## extract TBF and g-prior info
    tbf <- attrs$distribution$tbf
    gPrior <- attrs$distribution$gPrior

    ## correct estimateMargLik option
    estimateMargLik <- estimateMargLik && (! tbf)

    ## correct MCMC option
    if(tbf)
    {
        mcmc <- McmcOptions(burnin=0L,
                            step=1L,
                            samples=sampleSize(mcmc))
    }
    
    ## compute the selected approximation to the marginal z posterior,
    ## including the log density and the random number generator for it
    marginalz <-
        if(useFixedZ)
        {
            ## this is for the special case with a fixed z
            list(logDens=function(z) ifelse(z == fixedZ, 0, -Inf),
                 gen=function(n=1) rep.int(fixedZ, n))
        } else {
            if(tbf && is(gPrior, "IncInvGammaGPrior"))
            {
                ## analytical solution possible in this case!
                ## compute posterior parameters of IncIG:
                a <- gPrior@a + (modelDim - 1) / 2               
                b <- gPrior@b + info$residualDeviance / 2
                
                list(logDens=
                     function(z){
                         logNormConst <-
                             ifelse(b > 0,
                                    a * log(b) - pgamma(b, a, log.p=TRUE) - lgamma(a),
                                    log(a))                         
                         logNormConst - (a + 1) * log1p(exp(z)) -
                             b / (1 + exp(z)) + z
                     },
                     gen=
                     function(n=1){
                         p <- runif(n=n)
                         
                         ret <-
                             if(b > 0)
                             {
                                 b / qgamma(p=(1 - p) * pgamma(b, a),
                                            shape=a) - 1
                             } else {
                                 (1 - p)^(-1/a) - 1
                             }
                         return(ret)
                     })
            } else {
                ## the usual way:
                getMarginalZ(info,
                             method=marginalZApprox,
                             verbose=verbose)
            }
        }
    
    ## pack the info in handy lists:
    
    ## pack the options
    options <- list(mcmc=mcmc,
                    estimateMargLik=estimateMargLik,
                    verbose=verbose,
                    debug=debug,
                    isNullModel=isNullModel,
                    useFixedZ=useFixedZ,
                    fixedZ=as.double(fixedZ),
                    useOpenMP=useOpenMP)

    ## start the progress bar (is continued in the C++ code)
    if(verbose)
    {
        cat ("0%", rep ("_", 100 - 6), "100%\n", sep = "")
    }
    
    ## then call C++ to do the rest:
    cppResults <- cpp_sampleGlm(
      model,
      attrs$data,
      attrs$fpInfos,
      attrs$ucInfos,
      attrs$fixInfos,
      attrs$distribution,
      attrs$searchConfig,
      options,
      marginalz)

    ## start return list
    results <- list()

    ## add info on tbf
    results$tbf <- tbf
    
    ## compute acceptance ratio
    results$acceptanceRatio <-
        if(tbf)
            1 ## we do Monte Carlo in the TBF case!
        else
            cppResults$nAccepted / mcmc@iterations

    ## are we again verbose?
    if(verbose)
    {
        if(tbf)
            cat("\nFinished MC simulation from approximate parameter posterior\n")
        else
            cat("\nFinished MCMC simulation with acceptance ratio",
                round(results$acceptanceRatio, 3),
                "\n")
    }

    ## compute log marginal likelihood estimate (which is only defined up to a constant, so
    ## it can only be used for model ranking (?))
    if(estimateMargLik)
    {
        results$logMargLik <-
            with(cppResults,
                 getLogMargLikEstimate(numeratorTerms=samples$margLikNumerator,
                                       denominatorTerms=samples$margLikDenominator, 
                                       highDensityPointLogUnPosterior=highDensityPointLogUnPosterior))

        ## append original things
        results$logMargLik <-
            c(results$logMargLik,
              list(numeratorTerms=cppResults$samples$margLikNumerator,
                   denominatorTerms=cppResults$samples$margLikDenominator, 
                   highDensityPointLogUnPosterior=cppResults$highDensityPointLogUnPosterior))         
    }


    ## start post-processing.

    ## abbreviation for the coefficients sample matrix (nCoefs x nSamples)
    simCoefs <-
        results$coefficients <- cppResults$samples$coefficients
     
    
    ## so the number of samples is:
    nSamples <- ncol(simCoefs)
    
    ## the linear predictor samples for the fitted data:
    fitted <- design %*% simCoefs

    # give some names
    
    rownames(results$coefficients) <- colnames(design)
    if(doGlm) rownames(results$coefficients)[1] <- "(Intercept)"
    
    ## samples from the predictive distribution for new data,
    ## if this is required
    predictions <-
        if(nNewObs)
        {
            ## get new data matrix
            tempX <- constructNewdataMatrix(object=object,
                                            newdata=newdata)
            
            ## copy old GlmBayesMfp object
            tempObj <- object
          
            ## correct model matrix in tempMod to new data matrix
            attr(tempObj, "data") <-
                list(x=tempX,
                     xCentered=scale(tempX, center=TRUE, scale=FALSE)) 
           
            if(correctedCenter==FALSE){
              ## so we can get the design matrix
              newDesign <- getDesignMatrix(modelConfig=config,
                                           object=tempObj,
                                           intercept=doGlm)
              
              ## so the linear predictor samples are
              linPredSamples <- newDesign %*% simCoefs
              
            } else if(correctedCenter==TRUE){
              
              newDesignUC <- getDesignMatrix(modelConfig=config,
                                             object=tempObj,
                                             intercept=doGlm,
                                             center=FALSE)
              
              oldDesignUC <- getDesignMatrix(modelConfig=config,
                                             object=object,
                                             intercept=doGlm,
                                             center=FALSE)
              
              oldMeans <- colMeans(oldDesignUC)
              
              start <- ifelse(doGlm, 2,1) #If there is an intercept we don't want to subtract mean
              for(k in start:length(oldMeans)){
                newDesignUC[,k] <- newDesignUC[,k] - oldMeans[k]
              }
              
              ## so the linear predictor samples are
              linPredSamples <- newDesignUC %*% simCoefs
            }
            
              linPredSamples
            
            
        } else {
            ## no prediction required.
            ## But we return an empty matrix because that is expected
            ## by the GlmBayesMfpSamples class!
            matrix(nrow=0L,
                   ncol=0L)
        }

    ## process all coefficients samples: fixed, FP, UC in this order.
    
    ## now we need the colnames of the design matrix:
    colNames <- colnames(attr(object,"data")$x)
    
    
    ## invariant: already coefCounter coefficients samples (rows of simCoefs) processed    
    coefCounter <- 0L

    ## the fixed covariate (and intercept samples)
    fixCoefs <- list ()
        
    if(doGlm){
      fixName <- attrs$termNames$fixed[1]
      mat <- simCoefs[1, ,
                      drop=FALSE]   
      coefCounter <- coefCounter + 1
      
      ## and also get the names
      rownames(mat) <- colNames[1]
      
      ## and this is located in the list
      fixCoefs[[fixName]] <- mat
    }
    
    ## start processing all fixed terms
    for (i in seq_along(fixList <- attrs$indices$fixed))
    {
      ## get the name of the fixed term
      fixName <- attrs$termNames$fixed[i+doGlm]
      
      ## check if this fixed covariate is included in the model
      if (i %in% config$fixTerms)
      {
        ## then we get the corresponding samples
        mat <- simCoefs[coefCounter +
                          seq_len (len <- length(fixList[[i]])), ,
                        drop=FALSE]   
        
        ## correct invariant
        coefCounter <- coefCounter + len
        
        ## and also get the names
        rownames(mat) <- colNames[fixList[[i]]]
        
        ## and this is located in the list
        fixCoefs[[fixName]] <- mat
      }
    }
            
   
    ## samples of fractional polynomial function values evaluated at grids will
    ## be elements of this list:
    bfpCurves <- list ()                  
    ## Note that only those bfp terms are included which also appear in the
    ## model configuration.
    
    ## start processing all FP terms
    for (i in seq_along (attrs$indices$bfp))
    {
        ## what is the name of this FP term?
        fpName <- attrs$termNames$bfp[i]

        ## the powers for this FP term
        p.i <- config$powers[[fpName]]
        
        ## if there is at least one power for this FP, we add a list element, else not.
        if ((len <- length(p.i)) > 0)
        {            
            ## determine additional grid values:
            obs <- attrs$data$x[, attrs$indices$bfp[i], drop = FALSE]

            ## if there is no grid in gridList, we take an additional scaled grid using the gridSize argument
            if (is.null(g <- gridList[[fpName]]))
            {
                g <- seq (from = min(obs),
                          to = max(obs),
                          length = gridSize)
            }
            
            ## the resulting total grid is:
            g <- union (obs, g)
            gridSizeTotal <- length (g)

            ## sort the grid
            g <- sort(g)

            ## and rearrange as column with name (needed for getFpTransforms)
            g <- matrix(g,
                        nrow = gridSizeTotal,
                        ncol = 1L,
                        dimnames = list (NULL, fpName))
            
            ## the part of the design matrix corresponding to the grid
            xMat <- getFpTransforms (g, p.i, center=TRUE)

            ## multiply that with the corresponding coefficients to get the FP curve samples
            mat <- xMat %*% simCoefs[coefCounter + seq_len (len), , drop = FALSE]

            ## correct invariant
            coefCounter <- coefCounter + len
            
            ## save the grid as an attribute of the samples
            attr(mat, "scaledGrid") <- g
            
            ## save position of observed values
            attr(mat, "whereObsVals") <- match(obs, g) 

            ## then write into list
            bfpCurves[[fpName]] <- mat
        }
    }

    ## uncertain fixed form covariates coefficients samples
    ## will be elements of this list:
    ucCoefs <- list ()
    ## Note that only those UC terms are included which also appear in the model
    ## configuration.
    

    ## start processing all UC terms
    for (i in seq_along(ucList <- attrs$indices$ucList))
    {
        ## get the name of the UC term
        ucName <- attrs$termNames$uc[i]

        ## check if this UC is included in the model
        if (i %in% config$ucTerms)
        {
            ## then we get the corresponding samples
            mat <- simCoefs[coefCounter +
                            seq_len (len <- length(ucList[[i]])), ,
                            drop=FALSE]   

            ## correct invariant
            coefCounter <- coefCounter + len

            ## and also get the names
            rownames(mat) <- colNames[ucList[[i]]]

            ## and this is located in the list
            ucCoefs[[ucName]] <- mat
        }
    }

    ## collect processed samples into S4 class object
    results$samples <- new("GlmBayesMfpSamples",
                           fitted=fitted,
                           predictions=predictions,
                           fixCoefs=fixCoefs,
                           z=cppResults$samples$z,
                           bfpCurves=bfpCurves,
                           ucCoefs=ucCoefs,
                           shiftScaleMax=attrs$shiftScaleMax,
                           nSamples=nSamples)
    
    ## finished post-processing.

    ## finally return the whole stuff.
    return(results)
}
