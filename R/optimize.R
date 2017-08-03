#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian fractional polynomials for GLMs
##        
## Time-stamp: <[optimize.R] by DSB Die 25/05/2010 16:20 (CEST)>
##
## Description:
## R code and interface to C++ code that could be used for the marginal likelihood
## approximation, for regression testing purposes.
##
## History:
## 08/12/2009   file creation
## 09/12/2009   add a wrapper which uses both cppOptimize and cppScalarBfgs
##              to reduce the number of iterations and get a good hessian at the
##              minimum
## 15/02/2010   Note that optimBfgs does not work well!
## 25/05/2010   remove the function which uses both optimize and bfgs;
##              move the bfgs interface to this file and rename it to "cppBfgs";
##              add roxygen doc
#####################################################################################


##' Interface to the internal C++ optimization routine "bfgs"
##'
##' @param x0 the start value
##' @param f_ the target function
##' @param min.x minimum bound on x (default \code{-Inf})
##' @param max.x maximum bound on x (default \code{+Inf})
##' @param prec precision (default \code{1e-5})
##' @param verbose be verbose? (not default)
##' @return A list with the following elements:
##' \describe{
##' \item{par}{the minimum abscissa found by the algorithm}
##' \item{inv.hessian}{the inverse Hessian at \code{par}}
##' \item{evaluations}{list of the function evaluation pairs: \code{args} and
##' \code{vals}} 
##' \item{code}{the convergence code. 0 is \dQuote{OK}, -1 is \dQuote{lost
##' precision}, and +1 is \dQuote{change not large enough}}
##' }
##' 
##' @keywords utilities internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
cppBfgs <- function(x0,                    
                          f_,                    
                          min.x=c(-Inf),         
                          max.x=c(+Inf),         
                          prec=1e-5,             
                          verbose=FALSE)         
{
    return(cpp_bfgs( as.double(x0),
                     f_,
                     as.double(min.x),
                     as.double(max.x),
                     as.double(prec),
                     verbose))
}



##' Interface to the internal C++ optimization routine "optimize"
##'
##' @param f_ the function
##' @param min.x  minimum bound on x
##' @param max.x maximum bound on x
##' @param prec precision (same default as for the original optimize function)
##' @return A list with the following elements:
##' \describe{
##' \item{par}{the minimum abscissa found by the algorithm}
##' \item{inv.hessian}{the inverse Hessian at \code{par}}
##' \item{evaluations}{list of the function evaluation pairs: \code{args} and
##' \code{vals}} 
##' }
##' 
##' @keywords utilities internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
cppOptimize <- function(f_,          
                        min.x,       
                        max.x,       
                        prec=.Machine$double.eps^0.25) 
{
  return(cpp_optimize(f_,
                      as.double(min.x),
                      as.double(max.x),
                      as.double(prec)))
}


