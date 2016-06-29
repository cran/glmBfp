#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[getMarginalZ.R] by DSB Mon 26/08/2013 17:19 (CEST)>
##
## Description:
## Construct a (smooth) marginal z density approximation from a model
## information list 
##
## History:
## 17/02/2010   split off from sampleGlm and also export the function
## 25/02/2010   Different strategy: now potentially all four methods are tried
##              for the approximation, and the user may specify the order.
##              The normal approximation is no longer a special case which issued
##              a warning; not it is just one approximation method among all.
##              First try "spline", then "logspline", because sometimes this can remedy
##              the "too close to 0" error with "spline".
## 12/04/2010   limit the computation time for pinv.new to prepare against
##              infinite loops in this function.
## 15/05/2010   implement a very robust density approximation, which is just the
##              linear interpolation of the saved points. This "linear" method
##              differs from the other methods in that it does not depend on the
##              Runuran generators.
## 17/05/2010   "linear" is now the default, because it is very robust and works
##              nicely!
## 25/05/2010   now the "logDensVals" may contain NaN's, which we do not want to
##              include in the construction of the approximation, so we discard
##              these pairs early enough.
#####################################################################################


##' Construct a (smooth) marginal z density approximation from a model
##' information list 
##'
##' @param info the model information list
##' @param method method for approximating the marginal density:
##' \describe{
##' \item{linear}{Linearly interpolate the points.}
##' \item{spline}{The saved points of the unnormalized density approximation are joined
##' by a \dQuote{monotonic} spline. The density is smoothed out to zero at the
##' tails. Since the spline might be slightly negative for extreme values, the positive part
##' is returned.}
##' \item{logspline}{The saved points of the log unnormalized density approximation are joined 
##' by a \dQuote{monotonic} spline, which is then exponentiated.}
##' \item{normalspline}{A \dQuote{monotonic} spline is fitted to the differences of the saved
##' log density values and the log normal approximation. The resulting spline function is
##' exponentiated and then multiplied with the normal density.}
##' \item{normal}{Just take the normal approximation.}
##' This may also be a vector with more than one method names, to select the modify the preference
##' sequence: If the first method does not work, the second is tried and so on. The normal
##' approximation always \dQuote{works} (but may give bad results).
##' }
##' @param verbose Echo the chosen method? (not default)
##' @param plot produce plots of the different approximation steps? (not default)
##' @return a list with the log of the normalized density approximation (\dQuote{logDens}) and 
##' the random number generator (\dQuote{gen}).
##'
##' @export
##' @keywords internal
getMarginalZ <- function(info,
                         method=
                         c("linear", "spline", "logspline",
                           "normalspline", "normal"),
                         verbose=FALSE,
                         plot=FALSE)
{   
    ## get the (ordered) z
    zVals <- info$negLogUnnormZDensities$args
    zOrder <- order(zVals)    
    zVals <- zVals[zOrder]

    ## the range of z
    zMin <- zVals[1L]
    zMax <- zVals[length(zVals)]
   
    ## get (roughly normalized) and ordered log density values
    logDensVals <- with(info,
                        - negLogUnnormZDensities$vals - logMargLik)
    logDensVals <- logDensVals[zOrder]

    ## remove pairs with NaN log density values
    isNaN <- is.nan(logDensVals)

    logDensVals <- logDensVals[! isNaN]
    zVals <- zVals[! isNaN]
    
    ## optionally plot that
    if(plot)
    {
        par(mfrow=c(2, 2))
        plot(zVals,
             exp(logDensVals),
             type="o",
             main="roughly normalized original values")
        abline(v=info$zMode)
    }
    
    ## the approximating normal density is:
    normalDens <- function(z, log=FALSE)
    {
        dnorm(x=z,
              mean=info$zMode,
              sd=sqrt(info$zVar),
              log=log)
    }

    ## show it if wished
    if(plot)
    {
        curve(normalDens,
              from=zMin, to=zMax,
              main="normal approximation",
              n=301L)
        abline(v=info$zMode)
    }


    ## decide on the range of the spline function
    constant <- 0.2 * (zMax - zMin)
    extendedZVals <- c(zMin - constant,
                       zVals,
                       zMax + constant)

    ## add zeroes to both ends of the density values
    ## to smooth out the approximation to zero.
    densVals <- c(0,
                  exp(logDensVals),
                  0)

    ## grid where the approximation is checked
    extendedZGrid <- seq(from=min(extendedZVals),
                         to=max(extendedZVals),
                         length=401L)

    ## what are the possible methods for computing the unnormalized z density?
    possibleMethods <-
        expression(
            ## --------------------
            spline=
        {
            ## now compute the spline
            cubicSplineFun <- splinefun(x=extendedZVals,
                                        y=densVals,
                                        method="monoH.FC")

            ## so the unnormalized z density approximation is
            function(z, log=FALSE)
            {
                ## set everything outside the range to zero
                ## (unfortunately this is not an option for "splinefun")
                ret <- ifelse(z > min(extendedZVals) & z < max(extendedZVals),
                              cubicSplineFun(z),
                              0)

                ## take the positive part as the spline may be slightly negative
                ## at the ends
                ret <- pmax(ret, 0)

                ## exponentiate them?
                if(log)
                    return(log(ret))
                else
                    return(ret)
            }
        },
            ## --------------------
            linear=
        {
            ## now compute the spline
            linearSplineFun <- approxfun(x=extendedZVals,
                                         y=densVals,
                                         method="linear",
                                         rule=2) # return 0 outside the range.

            ## so the unnormalized z density approximation is
            function(z, log=FALSE)
            {
                ret <- linearSplineFun(z)

                ## exponentiate them?
                if(log)
                    return(log(ret))
                else
                    return(ret)
            }
        },
            ## --------------------            
            logspline=
        {
            ## add very small values to both ends of the log density values
            ## to smooth out the spline
            logZero <- min(logDensVals) - 1e5
            extendedLogDensVals <- c(logZero,
                                     logDensVals,
                                     logZero)

            ## now compute the spline
            cubicSplineFun <- splinefun(x=extendedZVals,
                                        y=extendedLogDensVals,
                                        method="monoH.FC")

            ## so the unnormalized z density approximation is
            function(z, log=FALSE)
            {
                ## set everything outside the range to zero
                ## (unfortunately this is not an option for "splinefun")
                logRet <- ifelse(z > min(extendedZVals) & z < max(extendedZVals),
                                 cubicSplineFun(z),
                                 logZero)

                ## exponentiate them?
                if(log)
                    return(logRet)
                else
                    return(exp(logRet))
            }
        },
            ## --------------------            
            normalspline=
        {            
            ## now fit the cubic spline to the differences of the log density values.
            ## In order to get smooth decrease to zero, we add zeroes to the log difference vector
            ## at both ends
            logDiffs <- c(0,
                          logDensVals - normalDens(zVals, log=TRUE),
                          0)

            cubicSplineFun <- splinefun(x=extendedZVals,
                                        y=logDiffs,
                                        method="monoH.FC")

            ## optionally plot the spline of the log differences
            if(plot)
            {
                curve(cubicSplineFun,
                      from=min(extendedZVals), to=max(extendedZVals),
                      n=301L,
                      main="spline of log differences")
                points(x=extendedZVals,
                       y=logDiffs)
                abline(v=info$zMode)
            }
            
            ## so the unnormalized z density approximation is
            function(z, log=FALSE)
            {
                ## set everything outside the range to zero
                ## (unfortunately this is not an option for "splinefun")
                splineVals <- ifelse(z > min(extendedZVals) & z < max(extendedZVals),
                                     cubicSplineFun(z),
                                     0)
                
                ## the log values:
                logRet <- splineVals + normalDens(z,
                                                  log=TRUE)

                ## exponentiate them?
                if(log)
                    return(logRet)
                else
                    return(exp(logRet))
            }
        },

            ## --------------------            
            normal=
        {
            normalDens
        })

    ## decide which method to try first  
    firstMethods <- match.arg(method, several.ok=TRUE)

    ## so what is the sequence this time?
    methodsSequence <- c(firstMethods,
                         setdiff(c("linear", "spline", "logspline", "normalspline", "normal"),
                                 firstMethods))
    
    ## then try the possible methods in this order
    for(m in methodsSequence)
    {
        ## get the unnormalized pdf
        unnormZdens <- eval(possibleMethods[[m]])

        ## first check that no infinite values occur!
        if(any(! is.finite(unnormZdens(extendedZGrid))))
        {
            warning("Infinite values for method ", m,
                    " in density approximation",
                    "going to the next method")
            
            ## if any is not finite, then go on to the next method
            next
        }

        ## then try to get the generator for that pdf
        generator <- try(getGenerator(method=m,
                                      unnormDensFun=unnormZdens,
                                      xVals=extendedZVals,
                                      unnormDensVals=densVals))
        
        ## if this resulted in an error (computational or time-out),
        ## go on to the next method,
        ## else break
        if(! inherits(generator, "try-error"))
        {
            if(verbose)
            {
                cat("\nTaking the", m, "approximation method\n")
            }
            break
        }
    }

    ## so the normalized log density is:
    logNormZdens <- function(z)            
    {
        unnormZdens(z, log=TRUE) - log(generator$normConst)
    }
    
    ## optionally plot the final normalized approximation
    if(plot)
    {    
        plot(extendedZGrid,
             exp(logNormZdens(extendedZGrid)),
             type="l",
             main="final approximation")
        abline(v=info$zMode)
    }

    ## return the list with the log density and the corresponding random number generator
    return(list(logDens=logNormZdens,
                gen=generator$generator))
}


##' Internal helper function which gets the generator (and normalizing
##' constant)
##'
##' @param method the methods string 
##' @param unnormDensFun the unnormalized density function
##' @param xVals the x values 
##' @param unnormDensVals the unnormalized density values
##' @return a list with the elements \dQuote{generator} and \dQuote{normConst},
##' containing the generator function (with argument \code{n} for the number
##' of samples) and the normalizing constant of the density, respectively.
##' 
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
getGenerator <- function(method,       
                         unnormDensFun,
                         xVals,
                         unnormDensVals)
{
    if(method=="linear")
    {
        ## we can build the generator ourselves.

        ## handy abbreviations
        nPoints <- length(xVals)
        z <- xVals
        y <- unnormDensVals
        
        ## first calculate the cdf points
        cdf <- numeric(nPoints)
        cdf[1L] <- 0
        for(i in 2:nPoints)
        {
            cdf[i] <- cdf[i - 1] +
                (y[i] + y[i - 1]) / 2 * (z[i] - z[i - 1])
        }

        ## normalize
        normConst <- cdf[nPoints]
        cdf <- cdf / normConst
        y <- y / normConst

        ## and now the corresponding inverse cdf (quantile function)
        quantFun <- function(p)
        {
            ## check argument
            stopifnot(0 <= p,
                      1 >= p)

            ## check the boundary cases
            if(p == 0)
                return(z[1])
            if(p == 1)
                return(z[nPoints])

            ## so now we are inside.

            ## find the right interval
            i <- 2L                             # we can start at 2
            while(p > cdf[i])
            {
                i <- i + 1L
            }
            ## so now cdf[i - 1] < p <= cdf[i]
            
            slope <- (y[i] - y[i - 1L]) / (z[i] - z[i - 1L])
            intercept <- y[i - 1L]
            const <- cdf[i - 1] - p

            solution <- (sqrt(intercept^2 - 2 * slope * const) - intercept) / slope
            
            return(solution + z[i - 1L])    
        }

        ## vectorize that
        quantFunVec <- Vectorize(quantFun)        
        
        ## so our generator function is
        gen <- function(n=1)
        {
            return(quantFunVec(runif(n=n)))
        }

        ## so we return this list:
        return(list(generator=gen,
                    normConst=normConst))        
        
    } else {
        ## some Runuran generator is built.

        ## the creation is somewhat slow,
        ## but the generation from it is very fast.
        
        ## However, sometimes pinv.new seems to get stuck in an
        ## infinite loop. Therefore we limit the computation time
        ## for this call (only).
        setTimeLimit(cpu=40)            # set time limit
        on.exit(setTimeLimit(cpu=Inf))  # delete time limit after exiting this
                                        # function. (either after normal exit or
                                        # after an error)
        
        unuranObject <- pinv.new(pdf=unnormDensFun,
                                 lb=min(xVals),
                                 ub=max(xVals),
                                 center=xVals[which.max(unnormDensVals)])
              
        ## get the normalizing constant
        ## (this is only possible *before* the object is packed!)
        normConst <- unuran.details(unuranObject, show=FALSE, return.list=TRUE)$area.pdf
    
        ## "pack" the object so that it can be saved just as a normal R object.
        unuran.packed(unuranObject) <- TRUE
    
        ## this function returns random variates from the pdf, where n is the number of variates to be
        ## produced
        gen <- function(n=1)
        {
            return(ur(unuranObject, n=n))
        }
    
        ## so we return this list:
        return(list(generator=gen,
                    normConst=normConst))        
    }
}
