#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[computeModels.R] by DSB Don 08/04/2010 13:48 (CEST)>
##
## Description:
## Compute model information (including the log marginal likelihood)
## for a given list of model configurations and glmBayesMfp output.
##
## History:
## 16/02/2010   file creation
## 08/04/2010   attribute writing works as intended.
#####################################################################################

##' @include helpers.R
{}

##' Compute model information for a given list of model configurations and glmBayesMfp output.
##'
##' If we want to compute the marginal likelihood and information necessary for
##' generating posterior samples for new models not encountered in the model search
##' done by \code{\link{glmBayesMfp}}, this function can be used: Provide it with the
##' models \code{configurations} to be interpreted in the context of the \code{object}
##' of class \code{\link{GlmBayesMfp}}. The result is again of the latter class, but contains
##' only the new models (similarly as the whole model space would consist of these and an
##' exhaustive search would have been conducted).
##' 
##' @param configurations list of the model configurations 
##' @param object the \code{\link{GlmBayesMfp}} object 
##' @param verbose be verbose? (default: only for more than 100 configurations)
##' @param debug be even more verbose and echo debug-level information? (not by default)
##' @return The \code{\link{GlmBayesMfp}} object with the new models. This can directly
##' be used as input for \code{\link{sampleGlm}}.
##' 
##' @export 
##' @keywords models regression
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
computeModels <- function(configurations,
                          object,
                          verbose=length(configurations) > 100L,
                          debug=FALSE)
{
    ## checks
    stopifnot(is.list(configurations),
              inherits(object, "GlmBayesMfp"),
              is.bool(verbose),
              is.bool(debug))
    
    ## get all attributes of the
    attrs <- attributes(object)

    ## include the logicals
    attrs$options$verbose <- verbose
    attrs$options$debug <- debug
    
    ## update the searchConfig argument to include the model configuration etc.
    ## so that the C++ side knows what to do!

    ## we do an "exhaustive" approach - only evaluate one model.
    ## So the C++ function glmExhaustive has the responsibility to
    ## handle this case.
    attrs$searchConfig$doSampling <- FALSE

    ## and these are the model configurations
    attrs$searchConfig$modelConfigs <- configurations
    
    ## then go C++
    result <- cpp_glmBayesMfp(attrs$data,
                              attrs$fpInfos,
                              attrs$ucInfos,
                              attrs$fixInfos,
                              attrs$searchConfig,
                              attrs$distribution,
                              attrs$options)

    ## C++ attaches the following attributes:

    ## numVisited
    ## inclusionProbs
    ## logNormConst

    ## we again just modify some of the attributes of the old object
    attrs$numVisited <- attr(result, "numVisited")
    attrs$inclusionProbs[] <- attr(result, "inclusionProbs")
    attrs$logNormConst <- attr(result, "logNormConst")
    attrs$names <- seq_along(result)
    
    ## so we can save much paperwork:
    attributes(result) <- attrs
    
    ## finally return the results
    return(result)
}
