#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* gmx *.* net]
## Project: Bayesian FPs for GLMs
##        
## Time-stamp: <[McmcOptions-methods.R] by DSB Die 25/05/2010 18:01 (CEST)>
##
## Description:
## Functions/methods for the MCMC formal class.
##
## History:
## 10/12/2009   copy from the master package
#####################################################################################

##' Compute the number of samples for a given MCMC options triple
##'
##' @param mcmcOptions the \code{\linkS4class{McmcOptions}} object
##' @return the resulting sample size
##'
##' @keywords programming
##' @export
sampleSize <-
    function(mcmcOptions)
{
    stopifnot(is(mcmcOptions, "McmcOptions"))

    return(as.integer(ceiling((mcmcOptions@iterations - mcmcOptions@burnin) / mcmcOptions@step)))
}


##' Convert samples to mcmc objects
##'
##' @param samples samples matrix or vector
##' @param mcmcOptions the \code{\linkS4class{McmcOptions}} object which was chosen for the
##' production of \code{samples}
##' @return an S3 class \dQuote{mcmc} object
##'
##' @importFrom coda mcmc
##' @keywords programming
##' @export
convert2Mcmc <- function(samples,
                         mcmcOptions)
{
    stopifnot(is(mcmcOptions, "McmcOptions"))
    
    return(coda::mcmc(data=samples,
                      start=mcmcOptions@burnin + 1L,
                      end=mcmcOptions@iterations,
                      thin=mcmcOptions@step))
}
