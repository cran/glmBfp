#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: 
##        
## Time-stamp: <[glmBfp-package.R] by DSB Mon 03/12/2012 18:19 (CET)>
##
## Description:
## Package description.
##
## History:
## 18/11/2009   file creation: copy and modify from master package.
## 08/12/2009   add cpp_optimize to the dynlibs
## 05/01/2010   add unuran.details to Runuran imports (since version 0.12, Josef
##              Leydold has implemented this with a return list where the normalizing
##              constant can be found.)
## 16/05/2010   add cpp_evalZdensity
## 25/05/2010   remove cpp_openmptest
## 27/07/2010   include import of "new" from methods
## 30/11/2012   add cpp_coxfit
## 03/12/2012   add "Cox" into the title
#####################################################################################

##' Bayesian inference for fractional polynomial models from the GLM and Cox
##' family 
##'
##' @name glmBfp-package
##' @aliases glmBfp
##' @docType package
##' @title Bayesian inference for fractional polynomial models from the GLM
##' and Cox family
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
##' @useDynLib glmBfp
##' cpp_glmBayesMfp cpp_bfgs cpp_optimize cpp_sampleGlm cpp_evalZdensity
##' cpp_coxfit
##' @importFrom graphics plot hist
##' @importFrom methods setClass setOldClass setGeneric setMethod representation
##' signature prototype initialize new
##' @importFrom Runuran pinv.new ur "unuran.packed<-" unuran.details
##' @keywords package
{}
