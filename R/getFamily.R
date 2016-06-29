#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[getFamily.R] by DSB Mit 06/02/2013 10:54 (CET)>
##
## Description:
## Helper for glmBayesMfp which extracts an S3 family object with additional elements
## from an option.
##
## History:
## 12/02/2010   split from glmBayesMfp.R file
## 12/05/2010   add cauchit link, remove inverse link from the accepted links
##              because we cannot deal with that in the generalized g-prior
##              framework.
## 06/02/2013   check for missing family, since there is no longer a default
##              in glmBayesMfp
#####################################################################################


##' Helper function for glmBayesMfp: Extracts an S3 family object
##'
##' Extracts an S3 family object, which (at least) the usual elements
##' "family", "link", "linkfun", "linkinv", "variance", "mu.eta", "dev.resids",
##' plus the additional elements:
##' "phi" which includes just the dispersion parameter to be used,
##' and "simulate" which can generate random variates for a given linear predictor and weight vector
##' (this function of course also uses the "phi" value)
##' 
##' @param family the family argument passed to \code{\link{glmBayesMfp}}
##' @param phi the dispersion argument passed to \code{\link{glmBayesMfp}}
##'
##' @return The returned family object also includes a custom \sQuote{init} function, 
##' which takes the response vector (or matrix) (\sQuote{y}) and the corresponding weight vector
##' (\sQuote{weights}), processes them to response vector \sQuote{y} and possibly altered weights
##' \sQuote{weights}, and includes starting values \sQuote{mustart} for the IWLS algorithm.
##' For example, here the binomial special case of two-column-response matrix is treated exactly in
##' the same way as in \code{\link{glm}}.
##'
##' @keywords internal
getFamily <- function(family,          
                      phi)                
{
    if(missing(family))
    {
        stop(simpleError("Option family must be specified!"))
    }
    
    ## check that dispersion is positive etc.
    phi <- as.numeric(phi)
    stopifnot(identical(length(phi),
                        1L),
              phi > 0)
    
    ## then process the family argument:

    ## convert character to function
    if (is.character(family))
        family <- get(family, mode = "function")
    ## call function to get object
    if (is.function(family))
        family <- family()
    ## check that family has correct form
    stopifnot(inherits(family, "family"),
              all(c("family", "link", "linkfun", "linkinv", "variance", "mu.eta", "dev.resids") 
                  %in% names(family)))
    
    ## check here if family and link are supported.
    ## (we only want to use links which always map into the valid mu space)

    ## first check the family,
    ## and along the way, we force the dispersion parameter
    ## and also choose which links are supported for that family.
    if(identical(family$family, "binomial"))
    {
        okLinks <- c("logit", "cauchit", "probit", "cloglog")
        family$phi <- 1
        family$simulate <- function(eta, weights)
        {
            ## check if weights are non-negative integers
            stopifnot(all((weights %% 1L) == 0L),
                      all(weights >= 0L)) 

            ## then simulate one observation for each lin pred value
            rbinom(n=length(eta),
                   size=weights,
                   prob=family$linkinv(eta))            
        }
    }
    else if(identical(family$family, "gaussian"))
    {
        okLinks <- c("log", "identity")
        family$phi <- phi
        family$simulate <- function(eta, weights)
        {
            rnorm(n=length(eta),
                  mean=family$linkinv(eta),
                  sd=sqrt(phi / weights))
        }
    }
    else if(identical(family$family, "poisson"))
    {
        okLinks <- c("log")
        family$phi <- 1
        family$simulate <- function(eta, weights)
        {
            ## here we do not use the weights at all!
            rpois(n=length(eta),
                  lambda=family$linkinv(eta))
        }
    }    
    else
    {
        stop(simpleError(paste("family",
                               family$family,
                               "not supported")))
    }


    ## warn if phi has changed
    if(family$phi != phi)
        warning(simpleWarning(paste("dispersion was changed from",
                                    phi,
                                    "to",
                                    family$phi)))

    ## then check the link
    if(! (family$link %in% okLinks))
    {
        stop(simpleError(paste("link",
                               family$link,
                               "not supported for family",
                               family$family)))
    }

    ## now attach a nicer initialize function, instead of just the expression e.g. in
    ## binomial()$initialize
    family$init <- function(y,          # the response vector or matrix
                            weights,    # the weights vector
                            nobs=length(weights)) # number of observations
    {
        ## create dummy variables which are sometimes needed in the
        ## initializer expression
        etastart <- start <- mustart <- NULL
        
        ## evaluate the initializer expression, which modifies y and weights
        ## and computes mustart values
        eval(family$initialize)

        ## return these
        return(list(y=y,
                    weights=weights,
                    linPredStart=family$linkfun(mustart)))
        ## use start values on the linear predictor scale for ease and clarity in C++.
    }
    
    return(family)
}
