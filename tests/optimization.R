#####################################################################################
## Author: Daniel Sabanés Bové [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: BFPs for GLMs.
##        
## Time-stamp: <[optimization.R] by DSB Die 25/05/2010 16:11 (CEST)>
##
## Description:
## Test the different optimization routines.
##
## History:
## 08/12/2009   start header;
##              add the new optimize function results and the number of function
##              evaluations to the results matrix
## 09/12/2009   also add combination of optimize and bfgs, which seems not a good idea
## 25/05/2010   remove testing of additional R optimization routines (see the archive
##              for the code if you need it) 
#####################################################################################


library(glmBfp)


## define simple univariate functions
funs <- list(function(x) log((x - 3)^2+ 0.1),
             function(x) (x - 3)^2,
             function(x) (x - 2)^4,
             function(x) abs(x - 3))

## and constraints
constraints <- c(min=-20, max=20)

## and plot range
xgrid <- seq(from=-30, to=30, length=201L)

## start value
start <- 0

par(mfrow=n2mfrow(length(funs)))

## now process all functions
for(f in funs)
{
    plot(xgrid, f(xgrid), type="l")

    cat("-------------------",
        "Now processing",
        paste(deparse(f), collapse=""),
        "-------------------\n")

    ## we can also compare with classic R code
    result0 <- optim(start, f, method="L-BFGS-B", hessian=TRUE,
                     lower=constraints["min"],
                     upper=constraints["max"])
    ## print(result0)
      
    ## and then with C++ scalar bfgs code
    result1 <- glmBfp:::cppBfgs(start, f, verbose=TRUE,
                                min.x=constraints["min"],
                                max.x=constraints["max"])
    ## print(result1)

    ## now the new optimize code
    result2 <- glmBfp:::cppOptimize(f,
                                    min.x=constraints["min"],
                                    max.x=constraints["max"])
    ## print(result2)

    ## assemble overview matrix
    results <- cbind(par=
                     c(result0$par,
                       result1$par,
                       result2$par),
                     inv.hessian=
                     c(1 / result0$hessian,
                       result1$inv.hessian,
                       result2$inv.hessian),
                     nEvaluations=
                     c(NA,
                       length(result1$evaluations$args),
                       length(result2$evaluations$args)))
    rownames(results) <- c("optim", "cppBfgs", "cppOptimize")

    ## and print it
    print(results)
}
