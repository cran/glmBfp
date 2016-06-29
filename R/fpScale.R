#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[fpScale.R] by DSB Die 16/02/2010 13:12 (CET)>
##
## Description:
## Scale an FP term variable appropriately before model search using "fpScale".
##
## History:
## 04/07/2008   copy from thesis function collection.
## 18/11/2009   document with Roxygen
#####################################################################################

##' Shift and scale a covariate vector (if wished) to have positive and small numbers.
##'
##' This function is (almost exactly) copied from package \dQuote{mfp} to be consistent.  
##'
##' @param x the covariate vector
##' @param scaling should shifting and scaling be performed (default)?
##'
##' @return list of \sQuote{shift} and \sQuote{scale} parameters appropriate for this
##' covariate vector.
##'
##' @keywords internal
fpScale <- function (x, scaling = TRUE)
{
    scale <- 1
    shift <- 0

    ## check scaling argument
    stopifnot(is.bool(scaling))
    
    if (isTRUE(scaling)) {
        if (min(x) <= 0) {
            z <- diff(sort(x))
            shift <- min(z[z > 0]) - min(x)
            shift <- ceiling(shift * 10)/10
        }
        range <- mean(x + shift)
        scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
    }
    return(list(shift = shift, scale = scale))
}

