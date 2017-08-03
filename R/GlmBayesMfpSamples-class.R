#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* gmx *.* net]
## Project: Bayesian FPs for GLMs
##        
## Time-stamp: <[GlmBayesMfpSamples-class.R] by DSB Mon 26/08/2013 17:30 (CEST)>
##
## Description:
## Formal class for samples from a single GlmBayesMfp model or a model average.
##
## History:
## 25/05/2010   file creation
## 28/05/2010   - change contents of "predictions": now only the linear predictors,
##              not predictive observations shall be saved here!
##              - rename "response" to "fitted", which shall also contain the
##              *linear predictors* for the fitted data (and not the *means*).
#####################################################################################

##' Class for samples from a single GlmBayesMfp model or a model average
##'
##' @note Note that for a model average, \code{nSamples} will typically *not*
##' be the number of samples available for each bfp curve, e.g.
##'
##' \describe{ 
##' \item{fitted}{fit samples on the linear predictor scale (\code{nObs} x
##' \code{nSamples})}
##' \item{predictions}{samples from the predictive distribution for new data
##' on the linear predictor scale, if \code{newdata} was provided
##' (\code{nrow(newdata)} x \code{nSamples})}
##' \item{fixCoefs}{the intercept and fixed covariate samples}
##' \item{z}{the log covariance factor samples}
##' \item{bfpCurves}{samples of fractional polynomial function values evaluated at
##' grids, contains one list element for each FP. Each element is
##' a matrix with the layout \code{nGridPoints x nSamples}.}
##' \item{ucCoefs}{uncertain fixed form covariates coefficients samples,
##' contains one list element for each fixed form covariate group.
##' Each element is a matrix with the layout \code{nCoefs x nSamples}.}
##' \item{shiftScaleMax}{transformation parameters}
##' \item{nSamples}{number of samples}
##' }
##'
##' @name GlmBayesMfpSamples-class
##' @keywords classes internal
setClass(Class="GlmBayesMfpSamples",
         representation=
         representation(fitted="matrix",
                        predictions="matrix",
                        fixCoefs="list",
                        z="numeric",
                        bfpCurves="list",
                        ucCoefs="list",
                        shiftScaleMax="matrix",
                        nSamples="integer"))
 


