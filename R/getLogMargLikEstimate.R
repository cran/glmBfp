#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[getLogMargLikEstimate.R] by DSB Mon 05/07/2010 13:36 (CEST)>
##
## Description:
## Compute the Chib-Jeliazkov log marginal likelihood estimate from the MCMC output
##
## History:
## 18/02/2010   separate code file;
##              correct bug for computing Omega(0): head(index, -0L) always gives
##              a length zero vector instead of the whole "index"!!!
#####################################################################################



##' Compute the Chib-Jeliazkov log marginal likelihood estimate from the MCMC output
##'
##' @param numeratorTerms terms for the sum in the numerator of the posterior density ordinate
##' estimate at the high density point
##' @param denominatorTerms terms for the sum in the denominator of the posterior density ordinate 
##' estimate at the high density point
##' @param highDensityPointLogUnPosterior log unnormalized posterior (i.e. likelihood times prior)
##' of the high density point
##' @param endLag up to which lag should the covariance incorporate the uncertainty? (default: 40)
##' @return list with resulting \dQuote{estimate} and \dQuote{standardError}.
##'
##' @keywords internal
##' 
##' @importFrom utils head tail
getLogMargLikEstimate <- function(numeratorTerms,
                                  denominatorTerms,
                                  highDensityPointLogUnPosterior,
                                  endLag=40L)
{
    ## note that here it is assumed that there are equally many numerator and
    ## denominator numbers, therefore
    stopifnot(identical(length(numeratorTerms),
                        length(denominatorTerms)),
              is.double(numeratorTerms),
              is.double(denominatorTerms),
              is.double(highDensityPointLogUnPosterior))

    ## first combine each numerator term with corresponding denominator term
    h <- rbind(numeratorTerms,
               denominatorTerms)

    ## the estimate for that is
    hHat <- rowMeans(h)
    
    ## so the estimate for the log posterior density ordinate is:
    logPosteriorDensityEstimate <- log(hHat[1]) - log(hHat[2])

    ## then compute the standard error for this estimate:
    hDiffHat <- h - hHat

    B <- ncol(hDiffHat)
    index <- seq_len(B)
    
    ## function to compute the autocovariance estimate for a given lag
    getOmega <- function(lag)
    {        
        ## get index for head part of time series
        headIndex <-
            if(identical(lag, 0L))
                index
            else
                head(index, - lag)

        ## and for tail
        tailIndex <-
            if(identical(lag, 0L))
                index
            else
                tail(index, - lag)

        ## now build up result
        ret <- matrix(data=0,
                      nrow=2L,
                      ncol=2L)

        for(l in seq_len(B - lag))
        {
            ret <- ret + tcrossprod(hDiffHat[, tailIndex[l]],
                                    hDiffHat[, headIndex[l]]) / B
        }

        return(ret)
    }

    ## start with Omega_0
    covMat <- getOmega(0L) / B

    ## then add other terms
    for(s in seq_len(endLag))
    {
        thisOmega <- getOmega(s)
        covMat <- covMat + (1 - s / (endLag + 1)) * (thisOmega + t(thisOmega)) / B
    }

    ## now covMat is ready to get delta ruled:
    aHat <- c(1 / hHat[1], - 1 / hHat[2])
    se <- sqrt(t(aHat) %*% covMat %*% aHat)
        
    ## return a list with the log marg lik estimate and the standard error
    return(list(estimate=highDensityPointLogUnPosterior - logPosteriorDensityEstimate,
                standardError=se))
}
