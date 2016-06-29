#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs.
## 
## Time-stamp: <[constructNewdataMatrix.R] by DSB Mit 06/01/2010 13:05 (CET)>
##
## Description:
## Internal function to construct the covariates matrix for new data based on the formula
## and scaling info in an existing GlmBayesMfp object, for passing it to prediction functions.  
##
## History:
## 10/12/2009   copy and modify from bfp package
#####################################################################################

##' Construct the covariates matrix for new data based on an existing GlmBayesMfp object
##'
##' This is an internal function, which constructs a model matrix for new covariate
##' data, based on the formula and scaling info in an existing GlmBayesMfp object.
##' The matrix can then be passed to prediction functions.
##'
##' @param object a valid \code{\link{GlmBayesMfp}}-Object
##' @param newdata the new covariate data as a data.frame (with the same
##' covariate names as in the call to \code{\link{glmBayesMfp}})
##' @return The (uncentered!) model matrix with the FP columns shifted and scaled as for the
##' original data.
##'
##' @keywords internal utilities
constructNewdataMatrix <- function(object,  
                                   newdata) 
{
    ## check class
    stopifnot(inherits(object, "GlmBayesMfp"))
    
    ## extract the covariates' formula and the names of the FP terms
    covariatesOnlyFormula <- attr (object, "formula")[-2]
    bfpTermNames <- attr (object, "termNames")$bfp

    ## get model matrix with the design variables for the new data
    newTerms <- terms (covariatesOnlyFormula, data = newdata) 
    ret <- model.matrix (newTerms, data = newdata)

    ## extract the transformation parameters
    tr <- attr (object, "shiftScaleMax")

    ## and scale the new FP columns like for the old data
    for (bfp in bfpTermNames)           
    {
        ## transform
        ret[, bfp] <- (ret[, bfp] + tr[bfp, "shift"]) / tr[bfp, "scale"]

        ## and afterwards check range: all transformed numbers must be positive!
        if(any(ret[, bfp] <= 0))
            stop(simpleError(paste("New data for covariate", bfp, "is out of range.")))
    }

    ## return the correctly transformed model matrix
    return(ret)
}



