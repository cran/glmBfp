##' Construct the design matrix for a given bfp GLM model
##'
##' This is an internal function to construct the UNCENTERED design matrix for a given
##' bfp GLM model.
##'
##' @param modelConfig the model configuration list which must have elements
##' \dQuote{powers} and \dQuote{powers}. Defaults to the configuration of the first element of 
##' @param object the \code{GlmBayesMfp} object, which is needed because it contains
##' the covariates matrix and indices vector
##' @param fixedColumns return the fixed columns inside the matrix (default) or not?
##' @return The design matrix, where the non-fixed part is columnwise centered (that is, the
##' colmeans are zero). 
##'
##' @keywords internal utilities
getUncenteredDesignMatrix <- function (modelConfig=object[[1]]$configuration, 
                                       object, 
                                       fixedColumns=TRUE)
{
  ## checks
  stopifnot(is.bool(fixedColumns))
  
  ## extract covariates matrices from object
  data <- attr (object, "data")
  full <- data$x
  fullCentered <- data$xCentered
  
  ## also extract index vector
  inds <- attr (object, "indices")
  
  ## then get powers / ucs from the model configuration
  powers <- modelConfig$powers
  ucSet <- modelConfig$ucTerms
  
  ## now the real code starts:
  
  ## check that only one covariate (the intercept) is fixed
  nFix <- length (inds$fixed)         
  stopifnot(identical(nFix,
                      1L))
  
  ## see how many ucs and fps are included
  ucColInds <- inds$uc %in% ucSet
  nUc <- sum(ucColInds)
  
  nFp <- length(unlist(powers))
  
  ## so we know how much space to reserve for the return matrix
  nColumns <-
    if(fixedColumns)
      nFix + nUc + nFp
  else
    nUc + nFp
  
  ret <- matrix (nrow = nrow (full),
                 ncol = nColumns)
  retColnames <- character (ncol (ret))
  
  ## invariant: already col columns written
  col <- 0                            
  
  if(fixedColumns)
  {
    ## fixed columns
    new <- full[, inds$fixed, drop = FALSE]
    newInds <- col + seq_along (inds$fixed)
    
    ret[, newInds] <- new
    retColnames[newInds] <- colnames (new)
    
    col <- col + nFix
  }
  
  ## fp part
  for (i in seq_along (inds$bfp)){
    pi <- powers[[i]]
    if (len <- length (pi)) {       # if there is at least one power
      new <- getFpTransforms (full[, inds$bfp[i], drop = FALSE], pi, center=FALSE)
      newInds <- col + seq_along (pi)
      
      ret[, newInds] <- new
      retColnames[newInds] <- colnames (new)
      
      col <- col + len
    }
  }
  
  ## uc part
  if (length (ucSet)){
    new <- full[, ucColInds, drop = FALSE]
    newInds <- col + seq_len (nUc)
    
    ret[, newInds] <- new
    retColnames[newInds] <- colnames (new)
    
    col <- col + nUc
  }
  
  ## attach dimnames 
  rownames (ret) <- rownames (full)
  colnames (ret) <- retColnames
  
  ## and then return the design matrix
  return (ret)
}