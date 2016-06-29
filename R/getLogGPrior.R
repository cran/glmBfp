#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[getLogGPrior.R] by DSB Mit 17/02/2010 15:55 (CET)>
##
## Description:
## Helper for glmBayesMfp which returns the normalized log g prior density.
##
## History:
## 12/02/2010   split from glmBayesMfp.R file
## 17/02/2010   reduce the lower bound for the integration to the smallest
##              value greater than zero instead of the machine epsilon, and compare the
##              difference to 1 with machine epsilon instead of the absolute error
##              (this is necessary e.g. for IG(0.001, 0.001))
#####################################################################################

##' Helper function for glmBayesMfp: Returns the normalized log g prior density
##'
##' Returns the normalized log of the prior density on g as an R function.
##' It is checked that the result is valid, that is that exp(f) integrates to 1 over the
##' positive real line.
##' 
##' @param gPrior the gPrior argument passed to \code{\link{glmBayesMfp}}
##'
##' @return The R function.
##'
##' @keywords internal
getLogGPrior <- function(gPrior)       
{
    ## get g prior choice
    if(isTRUE(is.character(gPrior)))
    {
        ## which of the choices has been selected by the user?
        gPrior <- match.arg(gPrior,
                            choices=("hyperg"))
        
        ## hyper-g prior with a = 4
        gPrior <-
            if(gPrior == "hyperg")
            {
                function(g) -2 * log1p(g)
            }        
    }
    
    ## check that the exp of the function it is a valid density function
    if(isTRUE(is.function(gPrior)))
    {
        XMIN <- .Machine$double.xmin
        EPS <- sqrt(.Machine$double.eps)
        
        integrand <- function(g) exp(gPrior(g))
        
        integral <- integrate(f=integrand,
                              lower=XMIN,
                              upper=Inf)
        
        if(integral$message != "OK")
            stop(simpleError(integral$message))
        if(abs(integral$value - 1) > EPS)
            stop(simpleError("the g-prior density must be proper and normalized"))
    }
    else
    {
        stop(simpleError("the g-prior argument must either be a valid string or a (vectorized) density function"))   
    }

    ## return the function
    return(gPrior)    
}
