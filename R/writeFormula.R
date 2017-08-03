#####################################################################################
## Author: Isaac Gravestock [isaac *.* gravestock *a*t* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Description:
## Construct a Survival model formula based on a glmBfp object
##
## History:
## 14/07/2015 Copied from CoxTBFs project functions-funsel-2.bma.R
#####################################################################################

{}

##' Construct a survival formula based on a glmBfp object with censInd not null.
##'
##' This is an internal function to construct a survival formula based on a glmBfp object
##' with censInd not null.
##'
##' @param models.listpart the glmBfp object for which to construct the survival formula
##' @param time.var the name of the time variable as character string
##' @param status.var the name of the censoring indicator variable as character string
##' @return The formula object based on the glmBfp object. 
##'
##' @keywords internal utilities

writeFormula <- function(models.listpart, time.var, status.var){
  surv.part <- paste("survival::Surv(", time.var, ", ", status.var, ") ~ ",sep="")
  
  #extract linear terms
  fix.parts <- models.listpart[[1]]$configuration$fixTerms
  fix.names  <- attributes(models.listpart)$termNames$fix
  
  uc.parts <- models.listpart[[1]]$configuration$ucTerms
  uc.names <- attributes(models.listpart)$termNames$uc
  
  shift <- attributes(models.listpart)$shiftScaleMax[,1]
  scale <- attributes(models.listpart)$shiftScaleMax[,2]
  
  if(length(models.listpart[[1]]$configuration$powers)>0){  
    #extract powers for terms which should have transformations
    powers <- models.listpart[[1]]$configuration$powers
    bfp.trans <- vector("character")
    for(i in 1:length(powers)){
      this.term <- ""
      if(length(powers[[i]])>0){
        
        this.name <- names(powers[i])
        this.power <- paste("c(",paste(powers[[i]],collapse=",", sep=""),")",sep="")
        this.term <- paste("fpTrans(",this.name,",", this.power, sep="")
        
        if(scale[i]!=0){
          this.term <- paste(this.term,", ",scale[i],sep="")
        }
        
        if(shift[i]!=0){
          this.term <- paste(this.term,", ",shift[i],sep="")
        }
        
        this.term <- paste(this.term, ")",sep="")
        
      }
      if(this.term != "") bfp.trans <- c(bfp.trans, this.term)
      
    }
    v <- paste(c(uc.names[uc.parts], bfp.trans,  fix.names[fix.parts]), collapse=" + ")
    return(formula(paste(surv.part,v)))
  } 
  
  #when there are no transformed variables
  v <- paste(c(uc.names[uc.parts], fix.names[fix.parts]), collapse=" + ")
  return(formula(paste(surv.part,v)))
}
