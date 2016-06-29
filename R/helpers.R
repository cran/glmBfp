#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs
##        
## Time-stamp: <[helpers.R] by DSB Die 16/02/2010 13:10 (CET)>
##
## Description:
## Helper functions
##
## History:
## 16/02/2010   file creation
#####################################################################################

##' Predicate checking for a boolean option
##'
##' @param x the object being checked
##' @return Returns \code{TRUE} if \code{x} is a length one logical vector (i.e., a
##' scalar)
##' 
##' @keywords internal
##' @author Daniel Sabanes Bove \email{daniel.sabanesbove@@ifspm.uzh.ch}
is.bool <- function(x)
{
    return(identical(length(x), 1L) &&
           is.logical(x))
}
