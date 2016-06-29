#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* gmx *.* net]
## Project: Bayesian FPs for GLMs
##        
## Time-stamp: <[McmcOptions-class.R] by DSB Mit 15/12/2010 13:57 (CET)>
##
## Description:
## Encapsulate the three canonical MCMC options in a formal class.
##
## History:
## 10/12/2009   copy from the master package
## 25/05/2010   include the optional argument "samples" in the constructor
## 08/07/2010   coerce default iterations to integer in constructor
## 13/12/2010   "samples" must have default value in the ctr
#####################################################################################

##' Class for the three canonical MCMC options
##'
##' The slots are:
##' \describe{
##' \item{iterations}{number of MCMC iterations}
##' \item{burnin}{number of burn-in iterations which are not saved}
##' \item{step}{only every step-th iteration is saved after the burn-in}
##' }
##'
##' @name McmcOptions-class
##' @keywords classes internal
setClass(Class="McmcOptions",
         representation=
         representation(iterations="integer",
                        burnin="integer",
                        step="integer"),
         validity=function(object){
             if(object@burnin < 0L )
             {
                 return("Burn-in must be non-negative")
             }
             else if(object@burnin >= object@iterations)
             {
                 return("Burn-in must be smaller than iterations")
             }
             else if(object@step < 1)
             {
                 return("Step size must be at least 1")
             }
         })
 
##' Constructor for class McmcOptions
##'
##' Note that the argument \code{samples} is included for convenience only -
##' you can specify it instead of \code{iterations}.
##' 
##' @param iterations number of MCMC iterations (default: \code{110,000})
##' @param burnin number of burn-in iterations which are not saved (default:
##' \code{10,000}) 
##' @param step only every step-th iteration is saved after the burn-in
##' (default: \code{10})
##' @param samples number of resulting samples (by default \code{10,000} will result)
##' @return the freshly built object of class \code{\linkS4class{McmcOptions}}
##'
##' @keywords classes
##' @export
McmcOptions <-
    function(iterations=as.integer(burnin + (step * samples)),
             burnin=1e4L,
             step=10L,
             samples=1e4L)
{
    return(new(Class="McmcOptions",
               iterations=iterations,
               burnin=burnin,
               step=step))
}
