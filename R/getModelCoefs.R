#####################################################################################
## Author: Isaac Gravestock [isaac *.* gravestock *a*t* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Description:
## Estimate shrunken coefficients from GlmBayesMfp object for Cox model
##
## History:
## 14/07/2015 Copied from CoxTBFs project functions-funsel-2.bma.R
#####################################################################################

##' @include sampleBma.R
##' @include getDesignMatrix.R
{}

##' Estimate shrunken coefficients from GlmBayesMfp object for Cox model
##'
##' This is an internal function to estimate shrunken coefficients from GlmBayesMfp object
##' for Cox models. It calls \code{\link{sampleBma}} and then calculates coefficents based
##' on a linear fit.
##'
##' @param models.listpart the glmBfp object for which to construct the survival formula
##' @param ... additional arguments to pass to \code{\link{sampleBma}}
##' @param sep should coefficients be returned separately (default=FALSE)
##' @return A named vector of coefficients. 
##'
##' @keywords internal utilities

getModelCoefs <- function(model.listpart,mcmc, sep=FALSE){
  
  DM <- getDesignMatrix(object=model.listpart, intercept=FALSE)
  n.coefs <- dim(DM)[2]
  
  if(n.coefs==0) return(0)
  
  if(sep==FALSE) {mcmc.obj <-McmcOptions()} else {mcmc.obj <- McmcOptions(samples = 200)}
  
  if(methods::hasArg(mcmc)) mcmc.obj <- mcmc
  
  model.bma <- sampleBma(model.listpart, mcmc=mcmc.obj,verbose=FALSE)
  
  
  #  print("break!")
  all.coefs <- lapply(1:ncol(model.bma$samples@fitted), function(i){  
    if(sum(is.na(model.bma$samples@fitted[,i]))>0) {
      print("This iteration from sampleBMA didn't work, return NA!")
      return(rep(NA, length(model.bma$samples@fitted[,i])))
    }
    result <- .lm.fit(x=as.matrix(DM), y=model.bma$samples@fitted[,i])$coefficients
    return(result)
  })
  
  if(sep==FALSE){
    if(n.coefs>1){ 
      arr <- array(unlist(all.coefs), dim=c(length(all.coefs[[1]]),length(all.coefs)))
      ret.all.coefs <- rowMeans(arr, na.rm=TRUE)
      #     ret.all.coefs <- rowMeans(all.coefs,na.rm=TRUE)
    } else ret.all.coefs <- mean(unlist(all.coefs), na.rm=TRUE)
    
    names(ret.all.coefs) <- colnames(DM)
    
    return(ret.all.coefs)
    
  } else if(sep==TRUE){
    return(all.coefs)
  }
  
}