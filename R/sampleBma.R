#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[sampleBma.R] by DSB Die 04/12/2012 09:29 (CET)>
##
## Description:
## Produce posterior samples from the Bayesian model average over models
## returned by glmBayesMfp, using MCMC.  
##
## History:
## 21/05/2010   file creation.
## 24/05/2010   first implementation (without testing)
## 25/05/2010   Do not return modelFreqs;
##              return a list with element "samples" from the S4 class
##              "GlmBayesMfpSamples".
##              test the first implementation.
## 28/05/2010   rename "response" to "fitted", which shall contain the
##              *linear predictors* for the fitted data (and not the *means*).
## 26/07/2010   insert many drop=FALSE statements in matrix subsetting
## 03/08/2010   there are always z samples (even in the null model), so we do
##              need to catch special cases.
## 26/11/2012   modifications to accommodate the TBF methodology
## 04/12/2012   modifications to accommodate the Cox models
#####################################################################################

##' @include GlmBayesMfp-methods.R
##' @include helpers.R
##' @include McmcOptions-class.R
##' @include McmcOptions-methods.R
##' @include sampleGlm.R
{}

##' Produce posterior samples from a Bayesian model average over GLMs / Cox
##' models
##'
##' Based on the result list from \code{\link{glmBayesMfp}}, sample from the
##' Bayesian model average (BMA) over the models contained in this list.
##'
##' If TBF methodology is used (which is specified within the \code{glmBayesMfp}
##' object), then Monte Carlo (MC) sampling is used. If the fully Bayesian,
##' generalized hyper-g prior methodology is used, then the sampling is based on
##' MCMC. Therefore, instead of only specifying the required number of samples
##' and the model probabilities, one also needs to specify the burn-in length
##' and the thinning parameter, which will be applied to every model from which
##' at least one sample is included in the average. Alternatively, you can ask
##' for MCMC marginal likelihood estimates for all models in the list. Then at
##' least \code{nMargLikSamples} will be produced for each model, whether
##' included in the BMA sample or not.
##'
##' @param object valid \code{GlmBayesMfp} object containing the models over
##' which to average
##' @param mcmc MCMC options object with class \code{\linkS4class{McmcOptions}},
##' specifying the number of required BMA samples (via \code{sampleSize(mcmc)}),
##' and the burn-in and thinning parameters applied to each model (see above).
##' If TBF is used, each sample is accepted, and the number of samples is given
##' by \code{\link{sampleSize}}(\code{mcmc}).
##' @param postProbs vector of posterior probabilites (will be normalized within
##' the function) for the weighting of the models in \code{object} (defaults to
##' the normalized posterior probabilities)
##' @param nMargLikSamples If this is non-\code{NULL}, it specified the number
##' of samples used for the marginal likelihood estimate for each model (see
##' above).
##' @param verbose should information on computation progress be given?
##' (default)
##' @param \dots optional further arguments already available for sampling from
##' a single model: \code{gridList}, \code{gridSize}, \code{newdata},
##' \code{weights}, \code{marginalZApprox}, \code{debug}, \code{useOpenMP}.
##' See \code{\link{sampleGlm}} for the meanings.
##' 
##' @return The result is a list with the following elements:
##' \describe{
##' \item{modelData}{data frame containing the result from the
##' \code{as.data.frame} function, and in addition BMA probabilities,
##' BMA frequencies in the sample, acceptance ratios of the MCMC
##' runs and optionally marginal likelihood estimates / standard
##' errors.}
##' \item{samples}{an object of S4 class \code{\linkS4class{GlmBayesMfpSamples}} 
##' containing the samples from the BMA.} 
##' }
##' 
##' @keywords models regression
##' @export
sampleBma <-
    function(object,
             mcmc=McmcOptions(),
             postProbs=posteriors(object),
             nMargLikSamples=NULL,
             verbose=TRUE,
             ...)
{
    ## Checks
    ## **************************************************
    
    ## check the object
    if(! inherits(object, "GlmBayesMfp"))
        stop(simpleError("object must be of class GlmBayesMfp"))

    nModels <- length(object)
    if(! (nModels >= 1L))
        stop(simpleError("there has to be at least one model in the object"))

    ## get attributes
    attrs <- attributes(object)
    
    ## check whether tbf is used
    tbf <- attrs$distribution$tbf

    ## is this a GLM or a Cox model?
    doGlm <- attrs$distribution$doGlm
    
    ## shall we do marginal likelihood estimations?
    estimateMargLik <- (! tbf) && (! is.null(nMargLikSamples))
    if(estimateMargLik)
    {
        ## then check if everything is OK
        stopifnot(is.numeric(nMargLikSamples),
                  nMargLikSamples > 100)
        nMargLikSamples <- as.integer(nMargLikSamples)
    } 
        
    ## other checks
    stopifnot(is.bool(verbose),
              is(mcmc, "McmcOptions"),
              identical(nModels,
                        length(postProbs)),
              postProbs >= 0)

    ## correct MCMC option
    if(tbf)
    {
        mcmc <- McmcOptions(burnin=0L,
                            step=1L,
                            samples=sampleSize(mcmc))
    }

    ## Start samples containers
    ## **************************************************

    ## "matrices"
    fitted <- matrix(nrow=attrs$data$nObs,
                     ncol=0L)

    ## we cannot give here the right number of rows,
    ## so use NULL to which we can cbind anything
    predictions <- NULL
            
    ## vectors
    z <- numeric()
    
    ## lists
    bfpCurves <-
      fixCoefs <-
      ucCoefs <- list()
    
    ## Distribute samples to models
    ## ************************************************** 

    ## determine the sample size
    nSamples <- sampleSize(mcmc)
              
    ## determine the model names and the normalized weights
    objNames <- as.numeric (names (object))
    postProbs <- postProbs / sum (postProbs)

    ## be sure that at least one model has positive weight
    stopifnot(any(postProbs > 0))

    ## draw model names
    modelFreqs <-
        if(identical(length(objNames),
                     1L))               # if there is only one model ...
        {
            rep(objNames, nSamples)
        }
        else                            # else more than one model ...
        {
            sample (objNames,
                    size = nSamples,
                    replace = TRUE,
                    prob = postProbs)
        }
    modelFreqs <- table (modelFreqs)
    nams <- names (modelFreqs)

    ## save model summary
    modelData <- as.data.frame(object)
    modelData[, c ("bmaProb", "bmaFreq")] <- 0
    modelData[, "bmaProb"] <- postProbs
    modelData[nams, "bmaFreq"] <- modelFreqs / nSamples

    ## prepare for acceptance rates
    modelData[, "acceptanceRatio"] <- 0
    
    ## prepare for marginal likelihood estimates
    if(estimateMargLik)
    {
        modelData[, c("margLikEstimate", "margLikError")] <- 0
    }
    
    ## Start sampling
    ## **************************************************

    ## echo sampling start
    if (verbose)
        cat ("\nStarting sampling ...")

    ## now sample from each model as often as indicated by the modelFreqs
    ## and save samples from it

    ## process every model in object 
    for (j in seq_along (object))
    {
        if (verbose)
            cat ("\nNow at model ", j, "...")

        ## get this model
        thisModel <- object[j]
        modName <- names(thisModel)
        
        ## determine if this model's samples are in BMA samples
        inBma <- modName %in% nams

        ## determine the number of samples we need from this model
        thisSampleSize <-
            if(inBma)
                ## if it is in the BMA, be sure that we have at least
                ## nMargLikSamples samples (if this is NULL, it counts as 0) 
                max(modelFreqs[modName],
                    nMargLikSamples)
            else
                ## else we only sample if we want marginal likelihood estimates 
                max(nMargLikSamples,
                    0L)
        
        ## decide if we need samples from this model
        if(thisSampleSize > 0L)
        {
           
            ## Sample from this model
            ## **************************************************

            ## adapt Mcmc object to reflect the correct number of samples
            thisMcmc <- McmcOptions(samples=thisSampleSize,
                                    burnin=mcmc@burnin,
                                    step=mcmc@step)

            ## and run the sampleGlm function
            thisOut <- sampleGlm(object=thisModel,
                                 mcmc=thisMcmc,
                                 estimateMargLik=estimateMargLik,
                                 verbose=verbose, # pass the verbose option
                                        # from this function
                                 ...) # and be sure to pass all other
                                        # options

            
            ## Save results
            ## **************************************************

            ## which samples do we keep from this model?
            keepSamples <-
                if(inBma)
                    ## the last modelFreqs[modName] samples
                    seq(from=thisSampleSize - modelFreqs[modName] + 1L,
                        to=thisSampleSize)
                else
                    ## none
                    integer()
                    
            ## save acceptance ratio
            modelData[modName, "acceptanceRatio"] <-
                thisOut$acceptanceRatio
            
            ## save the marginal likelihood estimates, if required
            if(estimateMargLik)
            {
                modelData[modName, c("margLikEstimate", "margLikError")] <-
                    unlist(thisOut$logMargLik[c("estimate", "standardError")])
            }

            ## save the predictive samples, if there are any
            if(nrow(thisOut$samples@predictions) > 0L)
            {
                predictions <- cbind(predictions,
                                     thisOut$samples@predictions[, keepSamples, drop=FALSE])
            }
            
            ## save z samples
            z <- c(z,
                   thisOut$samples@z[keepSamples])
            
            ## save all fixed coefs samples
            for(fixName in names(thisOut$samples@fixCoefs))
            {
              ## get these samples
              thisFixSamples <- thisOut$samples@fixCoefs[[fixName]]
              
              ## append these samples
              fixCoefs[[fixName]] <-
                cbind(fixCoefs[[fixName]],
                      thisFixSamples[, keepSamples, drop=FALSE])
            }    
            
            
            

            ## save the fit samples on the linear predictor scale
            fitted <- cbind(fitted,
                            thisOut$samples@fitted[, keepSamples, drop=FALSE])

            ## save all bfp curve samples
            for(bfpName in names(thisOut$samples@bfpCurves))
            {
                ## get these samples
                thisBfpSamples <- thisOut$samples@bfpCurves[[bfpName]]

                ## append these samples
                bfpCurves[[bfpName]] <-
                    cbind(bfpCurves[[bfpName]],
                          thisBfpSamples[, keepSamples, drop=FALSE])

                ## copy the grid attributes
                ## (actually we would not need to copy them every time)
                attr(bfpCurves[[bfpName]], "scaledGrid") <-
                    attr(thisBfpSamples, "scaledGrid")
                attr(bfpCurves[[bfpName]], "whereObsVals") <-
                    attr(thisBfpSamples, "whereObsVals")

            }

            ## save all uc coefs samples
            for(ucName in names(thisOut$samples@ucCoefs))
            {
                ## get these samples
                thisUcSamples <- thisOut$samples@ucCoefs[[ucName]]

                ## append these samples
                ucCoefs[[ucName]] <-
                    cbind(ucCoefs[[ucName]],
                          thisUcSamples[, keepSamples, drop=FALSE])
            }            
        }
    }

    ## Collect things in return list
    ## **************************************************

    ## be sure that predictions is a matrix, even if we have no
    ## predictive samples
    predictions <-
        if(is.null(predictions))
            matrix(nrow=0, ncol=0)
        else
            predictions
        
    ret <- list(modelData=modelData,
                samples=
                new("GlmBayesMfpSamples",
                    fitted=fitted,
                    predictions=predictions, 
                    fixCoefs=fixCoefs,
                    z=z,
                    bfpCurves=bfpCurves,
                    ucCoefs=ucCoefs,
                    shiftScaleMax=attrs$shiftScaleMax,
                    nSamples=nSamples))

    
    ## finally return the whole stuff.
    return(ret)
}


