#####################################################################################
## Author: Isaac Gravestock [isaac *.* gravestock *a*t* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Description:
## Transform variables in formulas with fractional polynomials
##
## History:
## 14/07/2015 Copy from CoxTBFs project
#####################################################################################


##' Transform formula variables
##'
##' Simple function to apply the Box Tidwell transformation to a variables in a formula.
##' Variable is first shifted and scaled
##' NewVar = (Var+shift)/scale
##' then transformed and optionally centered.
##' Can be used in formulas as poly() is.
##'
##' @param var the variable to transform
##' @param powers one or more powers
##' @param scale value to scale the variable after shifting (default=1)
##' @param shift value to shift the variable by (default=0)
##' @param center center the variable when tranforming.
##' 
##' @return the transformed vector
##' @keywords utilities
##' @export
##' 

fpTrans <- function(var, powers=1,scale=1, shift=0, center=TRUE ){
  varname <- deparse(substitute(var))
  
  varm<- as.matrix(var)
  
  newname <- if(shift !=0) {paste("(",varname,"+", shift,")",sep="")} else varname
  if(scale!= 1) newname <- paste(newname,"/",scale,sep="")
  
  colnames(varm) <- newname
  
  if(scale==0) stop("Attempting to scale (divide) by 0 in fpTrans")
  
  varm <- (varm+shift)/scale;
  
  #   out <- matrix(nrow=length(var), ncol=length(powers))
  
  res <- getFpTransforms(vec = varm, powers = powers, center = center)
  return(res)
}