#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[getFpTransforms.R] by DSB Fre 15/06/2012 16:05 (CEST)>
##
## Description:
## Transform a variable according to FP transformation formula and attach proper
## names to the resulting design matrix.
##
## History:
## 06/01/2010   copy from package bfp, slightly pimp and add Roxygen stuff
## 15/06/2012   simplify %bt% documentation for new roxygen package
#####################################################################################


##' Box Tidwell transformation
##'
##' Simple function to apply the Box Tidwell transformation to a vector, when a (single) power is
##' provided.
##'
##' @param x the numeric vector
##' @param pow the single power
##' @return the transformed vector
##'
##' @name boxTidwell
##' @aliases boxTidwell '%bt%'
##' @keywords utilities internal
'%bt%' <- function (x,            
                    pow)
{
    ## checks
    stopifnot(identical(length(pow),
                        1L),
              is.numeric(x),
              is.numeric(pow))

    ## transformation
    if (pow != 0)
    {
        x^pow
    }
    else ## pow == 0
    {
        log (x)
    }
}


##' Get the FP transforms matrix of a given covariate vector
##'
##' Get the (centered) FP transforms matrix of a given covariate vector, when the corresponding 
##' power vector (with at least one element) is provided.
##'
##' @param vec positive (== already shifted and scaled) column vector with proper colname
##' @param powers power vector with at least one element
##' @param center center the columns of the FP transforms matrix around zero? (default)
##' @return the FP transforms matrix with proper colnames.
##'
##' @keywords utilities internal
getFpTransforms <- function (vec,      
                             powers,   
                             center=TRUE)
{
    ## internal helper function, only required here:
    getTransformName <- function (name, pow)
    {
        ## checks
        stopifnot(identical(length(pow),
                            1L),
                  is.character(name))

        ## create name
        if(pow != 0)
        {
            paste (name, pow,
                   sep = "^")
        }
        else ## pow == 0
        {
            paste ("log(", name, ")",
                   sep = "")
        }
    }

    ## get properties of arguments
    name <- colnames (vec)
    vec <- as.vector (vec)
    len <- length (vec)
    np <- length (powers)

    ## start return matrix
    ret <- matrix (nrow = len, ncol = np)
    retColnames <- character (np)

    ## begin of recursion:

    ## start with vector of ones
    lastCol <- 1

    ## and power "log"
    lastPow <- 0

    ## invariant: already numRepPow times was this power repeated
    numRepPow <- 0L                     

    ## invariant: about to write ith column
    for (i in seq_along(powers))
    {
        ## extract this power into "pi",
        ## check if it is the same as the last power and (afterwards) if we are *not* at the first power.
        if ((pi <- powers[i]) == lastPow && i != 1L)
        {
            ## repeated powers case
            numRepPow <- numRepPow + 1L

            ## so we apply the rule for repeated powers
            lastCol <- lastCol * log(vec)            

            ## the name is a bit complicated
            retColnames[i] <-           
                if (pi == 0)
                {
                    ## log was repeated
                    paste (getTransformName (name, 0),
                           numRepPow + 1L,
                           sep = "^")
                }
                else
                {
                    ## other power was repeated
                    tmp <- paste (getTransformName (name, pi),
                                  getTransformName (name, 0),
                                  sep =  "*")
                    
                    if (numRepPow > 1L)
                        tmp <- paste (tmp, numRepPow, sep = "^")
                    
                    tmp
                }
        }
        else
        {
            ## normal case without repeated power:
            
            ## vector results from Box Tidwell transform
            lastCol <- vec %bt% pi

            ## and the name is easy
            retColnames[i] <- getTransformName (name, pi)

            ## reset power invariants
            numRepPow <- 0L
            lastPow <- pi
        }

        ## write the last column
        ret[, i] <- lastCol
    }

    ## optionally center the transforms matrix
    if(isTRUE(center))
    {
        ret <- scale(ret, center=TRUE, scale=FALSE)
    }

    ## attach proper colnames
    colnames (ret) <- retColnames

    ## and finally return that.
    return (ret)
}



