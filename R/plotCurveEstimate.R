#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[plotCurveEstimate.R] by DSB Mit 23/01/2013 18:14 (CET)>
##
## Description:
## Plot predictor curve estimates for a given samples object, produced by sampleGlm.
##
## History:
## 07/01/2010   copy and modify from package bfp (no method, add Roxygen chunks etc)
## 25/05/2010   expect input of S4 class "GlmBayesMfpSamples";
##              change expected layout of samples matrices (now nParameters x
##              nSamples)
## 03/08/2010   rug must be painted after the matplot call
#####################################################################################

##' @include hpds.R
{}

##' Function for plotting a fractional polynomial curve estimate
##'
##' Plot a fractional polynomial curve estimate using samples from a single
##' GLM / Cox model or a model average. 
##'
##' @param samples an object of class \code{\linkS4class{GlmBayesMfpSamples}},
##' produced by \code{\link{sampleGlm}} and \code{\link{sampleBma}}.
##' @param termName string denoting an FP term, as written by the
##' \code{\link[=as.data.frame.GlmBayesMfp]{as.data.frame}} method
##' @param plevel credible level for pointwise HPD, and \code{NULL} means
##' no pointwise HPD (default: 0.95). The pointwise intervals are plotted in
##' blue color.
##' @param slevel credible level for simultaneous credible band (SCB),
##' \code{NULL} means no SCB (defaults to \code{plevel}). The simultaneous
##' intervals are plotted in green color.
##' @param plot if \code{FALSE}, only return values needed to produce the
##' plot, but do not plot (default is \code{TRUE}, so a plot is made)
##' @param rug add a rug to the plot? (default: \code{FALSE})
##' @param addZeros include zero samples for models where the covariate is not
##' included? (default: \code{FALSE}) If \code{TRUE}, this changes the
##' interpretation of the samples, and therefore curve estimates based on these
##' samples: it is no longer conditional on inclusion of the covariate, but
##' marginally over all models, also those not including the covariate.

##' @param \dots further arguments for plotting with \code{\link{matplot}}
##' @return a list of various plotting information:
##' \item{original}{grid on the original covariate scale}
##' \item{grid}{grid on the transformed scale}
##' \item{mean}{pointwise mean curve values}
##' \item{plower}{lower boundaries for pointwise HPD}
##' \item{pupper}{upper boundaries for pointwise HPD}
##' \item{slower}{lower boundaries for SCB}
##' \item{supper}{upper boundaries for SCB}
##' \item{obsVals}{observed values of the covariate on the original scale}
##' \item{partialResids}{not implemented: partial residuals}
##' \item{transform}{vector of shift and scale parameter}
##'
##' @keywords regression
##' 
##' @importFrom graphics abline curve matplot par
##' @importFrom methods is
##' 
##' @export
plotCurveEstimate <-
    function (samples,                    
              termName,                 
              plevel = 0.95,            
              slevel = plevel,          
              plot = TRUE,              
              rug=FALSE,
              addZeros=FALSE,
              ...)
{
    ## check class of "samples"
    stopifnot(is(samples, "GlmBayesMfpSamples"))

    ## extract samples and its attributes (they might be overwritten by cbind
    ## below) 
    mat <- samples@bfpCurves[[termName]]
    attrs <- attributes(mat)
    
    ## check that there are samples for this covariate
    if (is.null(mat))
        stop ("There were no samples which include ", termName, " in this model sample!\n")   
    ## add zeros? todo: this should also work if no samples include this covariate!
    if(addZeros)
    {
        mat <- cbind(mat,
                     matrix(data=0,
                            nrow=nrow(mat),
                            ncol=(samples@nSamples - ncol(mat))))
    }    
    
    ## start return list
    ret <- list ()

    ## x values
    ret$grid <- g <- as.vector (attrs$scaledGrid)
    tr <- samples@shiftScaleMax[termName, c ("shift", "scale")]
    ret$original <- g * tr[2] - tr[1]

    ## compute pwise data
    ret$mean <- rowMeans(mat, na.rm=TRUE)

    if (! is.null(plevel))
    {
        plowerUpper <- apply(mat, 1, empiricalHpd, level = plevel)
        ret$plower <- plowerUpper[1, ]
        ret$pupper <- plowerUpper[2, ]
    }

    ## simultaneous credible band around the mean
    if (! is.null(slevel))
    {
        bandData <- scrHpd(mat, level = slevel, mode = ret$mean)
        ret$slower <- bandData[, "lower"]
        ret$supper <- bandData[, "upper"]
    }

    ## todo: add generalized partial residuals?
    
    ## ## partial residuals, attention because of possible ties between observed grid values in data!
    ## resids <- residuals (model)
    
    pos <- attrs$whereObsVals

    ## partialResids <- ret$mean[pos] + resids
    partialResids <- NULL
    
    if (plot){
        ## determine plotting arguments for matlines
        matplotList <- list (...)
        if (is.null (matplotList$xlab))
            matplotList$xlab <- termName
        if (is.null (matplotList$ylab)){
            front <- paste ("Average partial predictor g(", termName, ")", sep = "")
            if (any (tr != c (0, 1))){
                middle <- " after the transform "
                back <-
                    if (tr[1] != 0){
                        if (tr[2] != 1)
                            paste(termName, "%<-% (", termName, " + ", tr[1], ") %/% ", tr[2])
                        else
                            paste(termName, "%<-%", termName, " + ", tr[1])
                    } else {
                        paste(termName, "%<-%", termName, "%/%", tr[2])
                    }
                annotation <- substitute (expression (paste (f, m, b)),
                                          list (f = front,
                                                m = middle,
                                                b = parse (text = back)[[1]]))
            } else {
                annotation <- front
            }
            matplotList$ylab <- eval (annotation)
        }
        if (is.null (matplotList$lty))
            matplotList$lty <- 1
        if (is.null (matplotList$col))
            matplotList$col <- c ("black", "blue", "blue", "green", "green")
        if (is.null (matplotList$type))
            matplotList$type <- "l"
        matplotList$x <- ret$original
        matplotList$y <- as.data.frame(ret[- match(c("original", "grid"), names(ret))])
        
        if (is.null (matplotList$ylim))
            matplotList$ylim <- range (c (partialResids, matplotList$y))

        ## and plot:

        ## first the points
        ret$obsVals <- ret$original[pos]
        ret$partialResids <- partialResids
        ## plot(ret$obsVals, ret$partialResids,
        ##      type="p",
        ##      xlab=matplotList$xlab,
        ##      ylab=matplotList$ylab,
        ##      ylim=matplotList$ylim,
        ##      cex = 0.5,
        ##      col = "gray")
       
        ## then the curves, so that they are not over painted over by points
        ## matplotList$add <- TRUE
        do.call (matplot, matplotList)

        ## possibly the rug
        rug <- as.logical(rug)
        if(isTRUE(rug))
        {
            rug (jitter (ret$obsVals), col = "gray")
        }

    }
    ret$transform <- tr

    invisible (ret)
}
