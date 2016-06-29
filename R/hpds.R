#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[hpds.R] by DSB Die 25/05/2010 21:32 (CEST)>
##
## Description:
## HPD functions
##
## History:
## 07/01/2010   copy + modify from package bfp
## 25/05/2010   change layout of samples argument in scrHpd
#####################################################################################

##' Construct an empirical HPD interval from samples
##'
##' Construct an empirical highest posterior density (HPD) interval from
##' samples which have been drawn from the distribution of a quantity of
##' interest.
##'
##' @param theta the vector of samples
##' @param level the credible level
##' @return  A vector with the estimated lower and upper bounds of the HPD
##' interval.
##'
##' @seealso \code{\link{scrHpd}}
##' @keywords htest
##' @export
empiricalHpd <- function (theta,
                          level)
{
    ## check that theta is numeric, and that level is in (0, 1)
    stopifnot(is.numeric(theta),
              0 < level && 1 > level)

    ## how many samples are saved in theta?
    nSamples <- length (theta)

    ## get the sorted samples vector
    thetaSorted <- sort.int (theta, method = "quick")

    ## how many different credible intervals with "level"
    ## do we need to compare?
    nIntervals <- ceiling (nSamples * (1 - level))

    ## these are the start indexes of the intervals
    startIndexes <- seq_len (nIntervals)

    ## and these are the end indexes of the intervals
    endIndexes <- nSamples - nIntervals + startIndexes
    
    ## which interval has the smallest range?
    smallestInterval <- which.min(thetaSorted[endIndexes] - thetaSorted[startIndexes])

    ## then return the bounds of this smallest interval
    return (c (lower=thetaSorted[startIndexes[smallestInterval]],
               upper=thetaSorted[endIndexes[smallestInterval]]))
}


##' Calculate an SCB from a samples matrix
##'
##' Calculate an SCB from a samples matrix, which minimizes
##' the absolute distances of the contained samples to a mode vector, at
##' each gridpoint. Therefore the SCB might be considered an \dQuote{HPD
##' SCB}.
##'
##' @param samples m by n matrix where m is the number of parameters,
##' n is the number of samples and  hence each (multivariate) sample is a column in
##' the matrix \code{samples}
##' @param mode mode vector of length m (defaults to the vector of medians)
##' @param level credible level for the SCB (default: 0.95)
##' @return A matrix with columns \dQuote{lower} and \dQuote{upper}, with the
##' lower and upper SCB bounds, respectively.
##'
##' @references Besag, J.; Green, P.; Higdon, D. \& Mengersen, K. (1995):
##'  \dQuote{Bayesian computation and stochastic systems (with
##'   discussion)}, \emph{Statistical Science}, 10, 3-66.
##' @seealso \code{\link{empiricalHpd}}
##' @keywords htest multivariate
##' @export
scrHpd <- function(samples,  
                   mode = apply (samples, 1, median), 
                   level = 0.95)
{
    ## extracts
    nPars <- nrow(samples)                  # the number of parameters
    nSamples <- ncol(samples)                  # the number of samples

    ## checks
    if (nPars != length (mode))
        stop ("mode vector must have same length as samples matrix!")
    stopifnot(level > 0 && level < 1)

    ## absolute distance from mode vector
    distance <- abs (sweep (samples, 1, mode)) 

    ## Calculate a simultaneous (k/ nSamples)*100% credible band
    ## using the ranks approach:    
    k <- floor (level * nSamples)

    ## rowwise (= parameterwise) ranks of distances
    rankdistance <- apply (distance, 1, rank) 

    ## maximum ranks in each multivariate sample
    tstari <- apply (rankdistance, 1, max)

    ## sort the maximum ranks
    ordtstari <- sort.int (tstari, method = "quick") 

    ## the required rank, which divides the samples in contained and rejected samples for the SCB,
    ## is: 
    tstar <- ordtstari[k]               

    ## now which vectors are inside the SCB?
    whichInside <- tstari <= tstar
    ## note that sum(whichInside) is possibly larger than k, so we have
    ## a larger empirical coverage of the resulting SCB.

    ## reduce the samples matrix accordingly.
    samples <- samples[, whichInside]

    ## the parameterwise ranges of these vectors form the SCB
    ## (just the convex hull of all sample vectors!) 
    ret <- t(apply(samples, 1, range))
    colnames(ret) <- c("lower", "upper")

    ## finally return the (nPars x 2) matrix
    return (ret)
}
