#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs.
## 
## Time-stamp: <[inclusionProbs.R] by DSB Die 01/12/2009 13:21 (CET)>
##
## Description:
## Compute (model averaged) posterior inclusion probabilites for the uncertain
## variables (including FP variables) based on a GlmBayesMfp object.
##
## History:
## 01/12/2009   modify from the bfp package, and add roxygen tags.
#####################################################################################

##' Compute posterior inclusion probabilites based on GlmBayesMfp object
##'
##' Compute (model averaged) posterior inclusion probabilites for the uncertain
##' variables (including FP variables) based on a GlmBayesMfp object.
##'
##' @param GlmBayesMfpObject the GlmBayesMfp object
##' @param postProbs the vector of posterior model probabilities, defaults to
##' the normalized probabilities in \code{GlmBayesMfpObject}
##' @return the resulting inclusion probabilities vector
##'
##' @export
##' @keywords utilities
inclusionProbs <- function (GlmBayesMfpObject,
                            postProbs=
                            posteriors(GlmBayesMfpObject,
                                       type="normalized"))
{
    postProbs <- postProbs / sum(postProbs)
    inds <- attr (GlmBayesMfpObject, "indices")
    termNames <- attr (GlmBayesMfpObject, "termNames")

    nams <- unlist (termNames[c ("bfp", "uc")])
    ret <- numeric (length (nams))
    names (ret) <- nams

    i <- 0

    for (j in seq_along (inds$bfp)){
        i <- i + 1
        present <- sapply (GlmBayesMfpObject, function (one) as.logical (length (one$configuration$powers[[j]])))
        ret[i] <- sum (present * postProbs)
    }

    for (j in seq_along (inds$ucList)){
        i <- i + 1
        present <- sapply (GlmBayesMfpObject, function (one) any (j == one$configuration$ucTerms))
        ret[i] <- sum (present * postProbs)
    }

    return (ret)
}

## ****************************************************************************************************
