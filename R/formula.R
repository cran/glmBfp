#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
##        
## Time-stamp: <[formula.R] by DSB Fre 15/06/2012 16:04 (CEST)>
##
## Description:
## Helper functions for formula construction.
##
## History:
## 18/11/2009   file creation: extra file for this seems to be better.
##              Try to get one Rd file with rdname roclets (see make.Rd2.roclet documentation 
##              in roxygen package)
## 15/06/2012   remove x documentation from uc so that it is not duplicate
#####################################################################################

##' Mark a covariate for transformation with fractional polynomials
##'
##' Using this function, you can mark a covariate for transformation with fractional polynomials
##' in the formula for \code{\link{glmBayesMfp}}.
##'
##' @param x the covariate name
##' @param max maximum degree for this FP
##' @param scale use pre-transformation shifting and scaling to avoid numerical problems?
##' @param rangeVals extra numbers if the shifting and scaling should consider values in this range
##'
##' @return name of the provided covariate, with the other function parameters as  attached attributes
##'
##' @keywords utilities
##' @export
##' @rdname formula
bfp <-function (x,                    
                max = 2,             
                scale = TRUE,         
                rangeVals = NULL)
{
    x <- deparse(substitute(x))
    attr(x, "max") <- max
    attr(x, "scale") <- scale
    attr(x, "rangeVals") <- rangeVals
    x
}

## ****************************************************************************************************

##' Mark a covariate expression for joint variable selection 
##'
##' Using this function, you can mark a covariate or a group of combined covariates for joint
##' variable selection (\dQuote{uncertain covariate fixed form covariate groups}) in the formula for
##' \code{\link{glmBayesMfp}}. 
##'
##' @return Just the name of the provided covariate
##'
##' @keywords utilities
##' @export
##' @rdname formula
`uc` <-
function (x)
{
    x <- deparse(substitute(x))
    x
}

