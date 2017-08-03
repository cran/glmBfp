#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[getDesignMatrix.R] by DSB Die 16/02/2010 13:13 (CET)>
##
## Description:
##  Construct the design matrix for a given model configuration and GlmBayesMfp object.
##
## History:
## 06/01/2010   copy from package bfp, slightly modify signature and add Roxygen comments 
## 07/07/2015   add center argument if we want the uncentered design matrix.
#####################################################################################

##' @include helpers.R
##' @include getFpTransforms.R
{}

##' Construct the design matrix for a given bfp GLM model
##'
##' This is an internal function to construct the (centered) design matrix for a given
##' bfp GLM model.
##'
##' @param modelConfig the model configuration list which must have elements
##' \dQuote{powers} and \dQuote{powers}. Defaults to the configuration of the first element of 
##' @param object the \code{GlmBayesMfp} object, which is needed because it contains
##' the covariates matrix and indices vector
##' @param intercept return the intercept column inside the matrix (default) or not?
##' @param center should the data be centered (default) or not?
##' @return The design matrix, where the non-fixed part is columnwise centered (that is, the
##' colmeans are zero). 
##'
##' @keywords internal utilities
getDesignMatrix <- function (modelConfig=object[[1]]$configuration, 
                             object, 
                             intercept=TRUE,
                             center = TRUE)
{
    ## checks
    stopifnot(is.bool(intercept))
    
    ## extract covariates matrices from object
    data <- attr (object, "data")
    full <- data$x
    fullCentered <- data$xCentered

    ## also extract index vector
    inds <- attr (object, "indices")

    ## then get powers / ucs from the model configuration
    powers <- modelConfig$powers
    ucSet <- modelConfig$ucTerms
    fixSet <- modelConfig$fixTerms
    
    ## now the real code starts:

    ## check that only one covariate (the intercept) is fixed
    nFix <- sum (inds$fixed!=0) 
    # stopifnot(identical(nFix, 1L))

    ## see how many ucs and fps are included
    ucColInds <- inds$uc %in% ucSet
    nUc <- sum(ucColInds)
    
    fixColInds <- inds$fix %in% fixSet
  

    nFp <- length(unlist(powers))

    ## so we know how much space to reserve for the return matrix
    nColumns <-
        if(intercept)
            1 + nFix + nUc + nFp
        else
            nFix + nUc + nFp
    
    ret <- matrix (nrow = nrow (full),
                   ncol = nColumns)
    retColnames <- character (ncol (ret))

    ## invariant: already col columns written
    col <- 0                            

    if(intercept)
    {
      # add the intercept
      ret[, 1] <- full[, 1, drop=FALSE]
      col <- col + 1
    }
    
    ## fp part
    for (i in seq_along (inds$bfp)){
        pi <- powers[[i]]
        if (len <- length (pi)) {       # if there is at least one power
            new <- getFpTransforms (full[, inds$bfp[i], drop = FALSE], pi, center=center)
            newInds <- col + seq_along (pi)

            ret[, newInds] <- new
            retColnames[newInds] <- colnames (new)

            col <- col + len
        }
    }

    ## uc part
    if (length (ucSet)){
        if(center) {
			new <- fullCentered[, ucColInds, drop = FALSE]
        } else if (!center) {
			new <- full[, ucColInds, drop = FALSE]
		}
        newInds <- col + seq_len (nUc)

        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)

        col <- col + nUc
    }

    
    # add any other fixed cols
    if(length(fixSet)){
      if(center){
        new <- fullCentered[, fixColInds, drop = FALSE]
      } else if (!center){
        new <- full[, fixColInds, drop = FALSE]
      }
      
      newInds <- col + seq_len (sum(fixColInds)) # nFix-1 because we already put the intercept in
      ret[, newInds] <- new
      retColnames[newInds] <- colnames (new)
      
      col <- col + nFix
    }
    
    ## attach dimnames 
    rownames (ret) <- rownames (full)
    colnames (ret) <- retColnames

    ## and then return the design matrix
    return (ret)
}



