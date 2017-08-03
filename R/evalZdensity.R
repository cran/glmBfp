#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[evalZdensity.R] by DSB Fre 29/07/2011 15:13 (CEST)>
##
## Description:
## Evaluate the unnormalized marginal z density in a given model at a number of
## z values. This is mainly an exploratory function to test new optimization /
## integration strategies via the R interface.
##
## History:
## 16/05/2010   file creation: modified from the starting point sampleGlm.R
## 18/05/2010   modify the interface: now we accept one single model configuration
##              list, and a GlmBayesMfp object in which it is interpreted.
##              This is better since we do not depend on an already fitted
##              model.
## 08/07/2010   add an option to get the *conditional* marginal density
##              f(y | z) instead of the usual f(y, z) which is interpreted
##              as an unnormalized f(z | y).
## 29/07/2010   add the new option to get a better Laplace approximation in the
##              case of binary logistic regression.
## 29/07/2011   now "higherOrderCorrection"
#####################################################################################

##' @include helpers.R
{}

##' Evaluate the (negative log) unnormalized marginal z density in a given
##' model. 
##'
##' Based on the result list from \code{\link{glmBayesMfp}}, for the first model
##' in the list, the marginal z density can be evaluated.
##'
##' @param config the configuration of a single \code{GlmBayesMfp} model. The
##' null model is not allowed. It is interpreted in the context of
##' \code{object}.  
##' @param object the \code{\link{GlmBayesMfp}} object 
##' @param zValues the z values
##' @param conditional return the approximate *conditional* density f(y | z)?
##' (not default)
##' @param debug print debugging information? (not default)
##' @param higherOrderCorrection should a higher-order correction of the
##' Laplace approximation be used? (not default)
##' 
##' @return the negative log marginal unnormalized density values at the
##' \code{zValues}. (Note the words \dQuote{negative}, \dQuote{log}, and
##' \dQuote{unnormalized} !!!) 
##' 
##' @keywords internal
##' @export
evalZdensity <-
    function(config,
             object,
             zValues,
             conditional=FALSE,
             debug=FALSE,
             higherOrderCorrection=FALSE)
{
    ## check the object
    if(! inherits(object, "GlmBayesMfp"))
        stop(simpleError("object must be of class GlmBayesMfp"))

    ## other checks
    stopifnot(is.list(config),
              is.numeric(zValues),
              is.bool(conditional),
              is.bool(debug),
              is.bool(higherOrderCorrection))
    
    ## get the old attributes of the object
    attrs <- attributes(object)

    ## check if the model is the null model
    modelSize <- length(unlist(config))
    isNullModel <- identical(modelSize, 0L)
    if(isNullModel)
        stop(simpleError("the null model has no z which could be sampled..."))
    
    ## pack the info in a handy list:
    
    ## pack the options
    options <- list(zValues=as.double(zValues),
                    conditional=conditional,
                    debug=debug,
                    higherOrderCorrection=higherOrderCorrection)

    ## then call C++ to do the rest:
    results <- cpp_evalZdensity(config,
                                attrs$data,
                                attrs$fpInfos,
                                attrs$ucInfos,
                                attrs$distribution,
                                options)
    
    ## finally return the whole stuff.
    return(results)
}
