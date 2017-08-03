#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[GlmBayesMfp-methods.R] by DSB Fre 07/12/2012 11:10 (CET)>
##
## Description:
## Additional convenience methods for GlmBayesMfp class objects.
##
## History:
## 12/02/2010   modify from bfp package
## 15/03/2010   adapt print.GlmBayesMfp to new S4 gPrior element
## 07/12/2012   change default for freq in as.data.frame.GlmBayesMfp
#####################################################################################

##' @include posteriors.R
{}

##' Extract method for GlmBayesMfp objects
##'
##' Extract a subset of models from a \code{\link{GlmBayesMfp}} object.
##'
##' 
##' @method [ GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param \dots transports the indexes of the models
##' @return The subsetted object.
##' @export 
##' 
##' @name Extract.GlmBayesMfp
##' @aliases Extract.GlmBayesMfp [.GlmBayesMfp
##' @keywords methods
##' @seealso \code{\link{glmBayesMfp}} 
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
'[.GlmBayesMfp' <- function (x, ...)       
{
    y <- NextMethod("[")
    mostattributes (y) <- attributes (x)
    names (y) <- names (x)[...]
    class(y) <- oldClass(x)
    y
}



##' Print a GlmBayesMfp object.
##'
##'  
##' @method print GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param \dots unused
##' @return Only used for its side effect
##' 
##' @export
##' 
##' @keywords methods
##' @seealso \code{\link{glmBayesMfp}} 
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
print.GlmBayesMfp <- function (x, ...)
{
    cat ("------------------------------ GlmBayesMfp-Output ------------------------------\n")
    cat (length (x), "multivariable fractional polynomial model(s) of total (visited/cached)",
         attr (x, "numVisited"), "for following covariates:\n\n")
    cat ("fixed:                ", paste(attr (x, "termNames")$fixed, collapse = ", "), "\n")
    cat ("uncertain fixed form: ", paste(attr (x, "termNames")$uc, collapse = ", "), "\n")
    cat ("fractional polynomial:", paste(attr (x, "termNames")$bfp, collapse = ", "), "\n")
    pr <- attr (x, "priorSpecs")
    cat ("\nPrior for g was,", class(pr$gPrior),
         "\nand a ", pr$modelPrior, " prior on the model space was used.",
         "\n")
}


##' Convert a GlmBayesMfp object into a data frame
##'
##' 
##' @method as.data.frame GlmBayesMfp
##' @param x valid \code{\link{GlmBayesMfp}} object
##' @param row.names optional rownames (default is to keep the names of the
##' \code{\link{GlmBayesMfp}} list) 
##' @param \dots unused
##' @param freq should empirical frequencies of the models in the sampling
##' path be given? (not default)
##' @return The data frame with the following columns:
##' \describe{
##'     \item{posterior}{the posterior model probabilities}
##'     \item{logMargLik}{the log marginal likelihood of the models}
##'     \item{logPrior}{the log prior probabilities of the models}
##' }
##' Additionally, for each uncertain fixed form covariates a column with the inclusion
##' status, and for each fractional polynomial a column with the powers are returned.
##' @seealso \code{\link{glmBayesMfp}} 
##' @keywords methods
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
##' 
##' @export 
as.data.frame.GlmBayesMfp <- function (x, 
                                       row.names = NULL,
                                       ...,        
                                       freq = FALSE
                                       )
{
    ## posterior probabilites:
    ret <- data.frame (posterior=
                       posteriors(x, type="normalized"))
    if (! is.null(attr(x, "searchConfig")$chainlength) &&
        freq)
        ret$frequency <- posteriors (x, type="sampling")

    ## row names
    row.names (ret) <-
        if(! is.null(row.names))
            row.names
        else
            names(x)    

    ## log marginal likelihood and log prior
    ret$logMargLik <- logMargLiks(x)
    ret$logPrior <- logPriors(x)

    ## fixed form covariates:
    for (i in seq_along (fixNames <- attr (x, "termNames")$fixed[-1]))
    {
      ret[[fixNames[i]]] <- sapply (x,
                                   function(one)
                                     i %in% one$configuration$fixTerms)
    }
    
    ## uncertain fixed form covariates:
    for (i in seq_along (ucNames <- attr (x, "termNames")$uc))
    {
        ret[[ucNames[i]]] <- sapply (x,
                                     function(one)
                                     i %in% one$configuration$ucTerms)
    }
    ## fractional polynomial:
    for (fpName in attr (x, "termNames")$bfp)
    {
        ret[[fpName]] <- sapply (x,
                                 function(one)
                                 paste(one$configuration$powers[[fpName]],
                                       collapse = ", "))
    }

    ret
}
