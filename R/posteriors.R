#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[posteriors.R] by DSB Die 25/05/2010 15:11 (CEST)>
##
## Description:
## Extract posterior model probability estimates from a GlmBayesMfp object.
##
## History:
## 01/12/2009   modify from bfp package, and add roxygen tags;
##              add logmarglik and logprior extractor functions
#####################################################################################

##' Extract posterior model probability estimates from a GlmBayesMfp object
##'
##' @param GlmBayesMfpObject the object
##' @param type type of posterior model probability estimates to be extracted from
##' \code{GlmBayesMfpObject}
##' @return the requested probs from all models
##'
##' @export
##' @keywords utilities
posteriors <- function (GlmBayesMfpObject,
                        type = c("normalized", "sampling"))
                        
{
    ## match the type argument
    type <- match.arg(type)

    ## index is translated from type
    ind <- switch(type,
                  normalized=1L,
                  sampling=2L)

    ## check ... 
    ## ... does the first model have the probability estimate?
    stopifnot(! is.na(GlmBayesMfpObject[[1]]$information$posterior[ind]))
    
    ## now extract the requested probs from all models
    sapply (GlmBayesMfpObject, function (one) one$information$posterior [ind])
}

## ****************************************************************************************************

##' Extract the log marginal likelihood estimates from a GlmBayesMfp object
##'
##' @param GlmBayesMfpObject the object
##' @return the vector of log marginal likelihood estimates
##'
##' @export
##' @keywords utilities
logMargLiks <- function(GlmBayesMfpObject)
{
    sapply(GlmBayesMfpObject,
           function(one) one$information$logMargLik)
}

## ****************************************************************************************************

##' Extract the log prior values from a GlmBayesMfp object
##'
##' @param GlmBayesMfpObject the object
##' @return the vector of log prior values
##'
##' @export
##' @keywords utilities
logPriors <- function(GlmBayesMfpObject)
{
    sapply(GlmBayesMfpObject,
           function(one) one$information$logPrior)
}




