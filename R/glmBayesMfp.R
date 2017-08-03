#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Time-stamp: <[glmBayesMfp.R] by DSB Mit 03/07/2013 22:57 (CEST)>
##
## Description:
## Main user interface for Bayesian inference for fractional polynomials in generalized linear
## models and Cox models. 
##
## History:
## 26/10/2009   file creation: copy and modify the old BayesMfp.R
## 27/10/2009   first only allow the binomial family, and write the code in an
##              extensible way, so that it later can be extended easily when binomial works
##              (11/12/2009 note: actually we now have generality for all GLM families.)
## 06/11/2009   include group sizes
## 09/11/2009   add Gauss Hermite quantiles and weights
## 18/11/2009   progress towards first testable version of the marginal likelihood
##              approximation,
##              use roxygen after cleaning the package from old hyper-g stuff.
## 01/12/2009   do not coerce model response to type double automatically,
##              because it could be a factor in the binomial case.
## 11/12/2009   rewrite passing of info to C++ and storage of attributes info,
##              to make it clearer and facilitate reuse in other functions, e.g. sampleGlm.
## 06/01/2010   extend getFamily to return self-crafted "simulate" function in the
##              return list. This is used e.g. by "sampleGlm".
## 12/02/2010   for clarity, split off helper functions into their separate files,
##              so we need to include them in the preamble
## 17/02/2010   split off null model information computation, and include the returned
##              list in the things passed to C++ (instead of only the log marginal
##              likelihood for the null model)
## 15/03/2010   also pass weights to C++ via the family list,
##              new g-prior class is expected in the prior list
## 12/04/2010   do not coerce totalNumber to integer but to double, as it is done in C++
##              in the same manner and keeps the number space large enough.
## 14/04/2010   add "useBfgs" and "largeVariance" options.
## 17/05/2010   useBfgs=FALSE ("optimize") is now the default because it is more
##              robust than Bfgs. Add "useOpenMP" option, which makes it easier
##              to switch than setting the environment variable OMP_NUM_THREADS.
## 21/05/2010   be more careful when passing "fpnames" to C++
## 25/05/2010   add "nObs" to "data" attribute of return list.
## 08/07/2010   add an empirical Bayes option which ranks the models in terms of
##              approximate *conditional* marginal likelihoods
## 29/07/2010   add the new option to get a better Laplace approximation in the
##              case of binary logistic regression.
## 08/07/2011   add the new modelPrior option "dependent"
## 29/07/2011   now "higherOrderCorrection"
## 15/02/2012   check that not two terms contain the same covariate in
##              the formula
## 21/11/2012   add "tbf" option
## 03/12/2012   add Cox regression
## 10/12/2012   fix bug: order according to survival times when Cox model is
##              requested, otherwise the computed deviances are wrong!!
## 06/02/2013   remove default for "family"
## 03/07/2013   - add offsets
##              - remove getNullModelInfo
## 26/05/2014   Added option (useFixedc) to calculate (or not) c factor using
##              mean of observations as in null model instead of alpha=0.
#####################################################################################

##' @include helpers.R
##' @include formula.R
##' @include fpScale.R
##' @include GPrior-classes.R
##' @include getFamily.R
{}


##' Bayesian model inference for fractional polynomial GLMs and Cox models 
##'
##' Bayesian model inference for fractional polynomial models from the generalized linear model
##' family or the Cox model is conducted by means of either exhaustive model space evaluation or posterior model
##' sampling. The approach is based on analytical marginal likelihood approximations, using
##' integrated Laplace approximation. Alternatively, test-based Bayes factors
##' (TBFs) are used.
##'
##' The formula is of the form \code{y ~ bfp (x1, max = 4) + uc (x2 + x3)}, that
##' is, the auxiliary functions \code{\link{bfp}} and \code{\link{uc}} must be
##' used for defining the fractional polynomial and uncertain fixed form
##' covariates terms, respectively. There must be an intercept, and no other
##' fixed covariates are allowed. All \code{max} arguments of the
##' \code{\link{bfp}} terms must be identical. \code{y} is the response vector
##' for GLMs or the vector of survival times for Cox regression. Note that Cox
##' regression is only implemented with TBFs.
##' 
##' The prior specifications are a list:
##' \describe{
##'   \item{gPrior}{A g-prior class object. Defaults to a hyper-g prior. See
##'   \code{\linkS4class{GPrior}} for more information.} 
##'   \item{modelPrior}{choose if a flat model prior (\code{"flat"}), a
##'   model prior favoring sparse models explicitly (default, \code{"sparse"}),
##'   or a dependent model prior (\code{"dependent"}) should be used.}
##' }
##' 
##' If \code{method = "ask"}, the user is prompted with the maximum
##' cardinality of the model space and can then decide whether to use
##' posterior sampling or the exhaustive model space evaluation.
##'
##' Note that if you specify only one FP term, the exhaustive model search
##' must be done, due to the structure of the model sampling algorithm.
##' However, in reality this will not be a problem as the model space will
##' typically be very small.
##' @param formula model formula
##' @param censInd censoring indicator. Default is \code{NULL}, but if
##' a non-\code{NULL} vector is supplied, this is assumed to be logical
##' (\code{TRUE} = observed, \code{FALSE} = censored) and Cox regression is
##' performed.  
##' @param data optional data.frame for model variables (defaults to the parent
##' frame)
##' @param weights optionally a vector of positive weights (if not provided, a
##' vector of one's)
##' @param offset this can be used to specify an _a priori_ known component to
##' be included in the linear predictor during fitting. This must be a numeric
##' vector of length equal to the number of cases (if not provided, a vector of
##' zeroes)
##' @param family distribution and link (as in the glm function). Needs to
##' be explicitly specified for all models except the Cox model.
##' @param phi value of the dispersion parameter (defaults to 1)
##' @param tbf Use TBF methodology to compute the marginal likelihood? (not
##' default) Must be \code{TRUE} if Cox regression is done.
##' @param empiricalBayes rank the models in terms of \emph{conditional}
##' marginal likelihood, using an empirical Bayes estimate of g? (not default)
##' Due to coding structure, the prior on g must be given in \code{priorSpecs}
##' although it does not have an effect when \code{empiricalBayes==TRUE}.
##' @param fixedg If this is a number, then it is taken as a fixed value of g,
##' and as with the \code{empiricalBayes} option, the models are ranked in terms
##' of conditional marginal likelihood. By default, this option is \code{NULL},
##' which means that g is estimated in a fully or empirical Bayesian way.
##' @param priorSpecs prior specifications, see details
##' @param method which method should be used to explore the  posterior model
##' space? (default: ask the user)
##' @param subset optional subset expression
##' @param na.action default is to skip rows with missing data, and no other
##' option supported at the moment 
##' @param verbose should information on computation progress be given?
##' (default)
##' @param debug print debugging information? (not default)
##' @param nModels how many best models should be saved? (default: 1\% of the
##' total number of (cached) models). Must not be larger than \code{nCache} if
##' \code{method == "sampling"}.
##' @param nCache maximum number of best models to be cached at the same time
##' during the model sampling, only has effect if method = sampling 
##' @param chainlength length of the model sampling chain (only has an effect if
##' sampling has been chosen as method) 
##' @param nGaussHermite number of quantiles used in Gauss Hermite quadrature
##' for marginal likelihood approximation (and later in the MCMC sampler for the
##' approximation of the marginal covariance factor density). If
##' \code{empiricalBayes} or a fixed g is used, this option has no effect.
##' @param useBfgs Shall the BFGS algorithm be used in the internal maximization
##' (not default)? Else, the default Brent optimize routine is used, which seems
##' to be more robust. If \code{empiricalBayes} or a fixed g is used, this
##' option has no effect and always the Brent optimize routine is used.
##' @param largeVariance When should the BFGS variance estimate be considered
##' \dQuote{large}, so that a reestimation of it is computed? (Only has an
##' effect if \code{useBfgs == TRUE}, default: 100)
##' @param useOpenMP shall OpenMP be used to accelerate the computations?
##' (default)
##' @param higherOrderCorrection should a higher-order correction of the
##' Laplace approximation be used, which works only for canonical GLMs? (not
##' default) 
##' @param fixedcfactor If TRUE sets the c factor assuming alpha is set to 0. Otherwise take alpha=mean(y)
##' @param  empiricalgPrior If TRUE uses the the observed isnformation matrix instead of X'X in the g prior. (Experimental)
##' @param centerX Center the data before fitting (FALSE)
##'
##' @aliases glmBayesMfp GlmBayesMfp
##' @return An object of S3 class \code{GlmBayesMfp}.
##' 
##' @keywords models regression
##' @export
glmBayesMfp <-
    function (formula = formula(data),
              censInd = NULL,
              data = parent.frame(),   
              weights,
              offset,
              family,       
              phi=1,
              tbf=FALSE,
              empiricalBayes=FALSE,
              fixedg=NULL,
              priorSpecs =            
              list(gPrior=HypergPrior(), 
                   modelPrior="sparse"),
              method = c ("ask", "exhaustive", "sampling"),
              subset,           
              na.action = na.omit,     
              verbose = TRUE,
              debug=FALSE,
              nModels,
              nCache=1e9,
              chainlength = 1e4,  
              nGaussHermite=20,
              useBfgs=FALSE,
              largeVariance=100,
              useOpenMP=TRUE,
              higherOrderCorrection=FALSE,
              fixedcfactor=FALSE,
              empiricalgPrior=FALSE,
              centerX=TRUE)
{
    ## checks
    stopifnot(is.bool(tbf),
              is.bool(verbose),
              is.bool(debug),
              is.bool(useBfgs),
              is.bool(empiricalBayes),
              is(priorSpecs$gPrior, "GPrior"),
              is.bool(useOpenMP),
              is.bool(higherOrderCorrection),
              is.bool(empiricalgPrior))

    ## see whether GLM or Cox is requested
    doGlm <- is.null(censInd)
    if(! doGlm)
    {
        ## and check censoring indicator vector in the Cox case.
        ## Also: we only have Cox regression with TBFs!
        stopifnot(is.logical(censInd),
                  tbf)

        ## set family to Gaussian to have some pseudo values in there
        family <- gaussian
    }

    ## see whether a fixed g is requested
    useFixedg <- ! is.null(fixedg)
    if(useFixedg)
    {
        ## then check whether it is a valid g
        ## and that there is no conflict with the empirical Bayes option 
        stopifnot(is.numeric(fixedg),
                  identical(length(fixedg), 1L),
                  fixedg > 0,
                  ! empiricalBayes)
    } else {
        fixedg <- 0
    }
    
    ## save call for return object
    call <- match.call()
    method <- match.arg (method)

    ## check and evaluate Gauss Hermite stuff
    nGaussHermite <- as.integer(nGaussHermite)
    gaussHermite <- statmod::gauss.quad(n=nGaussHermite, kind="hermite")
    
    ## evaluate family, this list then also includes the dispersion
    family <- getFamily(family, phi)    

    ## get model prior choice
    priorSpecs$modelPrior <- match.arg(priorSpecs$modelPrior,
                                       choices=c("flat", "sparse", "dependent"))
    
    ## evaluate call for model frame building
    m <- match.call(expand.dots = FALSE)

    ## select normal parts of the call
    temp <- c("", "formula", "data", "weights", "offset", "subset", "na.action") # "" is the function name
    m <- m[match(temp, names(m), nomatch = 0)]
   
    ## sort formula, so that bfp comes before uc
    ## filter special parts in formula: uncertain covariates (uc) and (Bayesian) fractional polynomials (bfp)
    special <- c("uc", "bfp")
        
    Terms <- if (missing(data))
            terms(formula, special)
        else
            terms(formula, special, data = data)

    tempVarNames <- rownames(attr(Terms, 'factors'))

    ## check if intercept is present
    if (! attr(Terms, "intercept"))
        stop(simpleError("there must be an intercept term in the model formula"))

    ucTempVarNames <- sort(tempVarNames[attr(Terms, "specials")$uc])
    bfpTempVarNames <- sort(tempVarNames[attr(Terms, "specials")$bfp])
    fixTempVarNames <- sort(tempVarNames[ - c(1,
                                              attr(Terms, "specials")$uc,
                                              attr(Terms, "specials")$bfp)])


    ## now sort the formula
    sortedFormula <- paste(deparse (Terms[[2]]),
                           "~ 1 +", 
                           paste(c(bfpTempVarNames, ucTempVarNames, fixTempVarNames),
                                 collapse = "+"))
    sortedFormula <- as.formula (sortedFormula)

    ## filter special parts in formula: uncertain covariates (uc) and (Bayesian) fractional polynomials (bfp)
    Terms <- if (missing(data))
        terms(sortedFormula, special)
    else
        terms(sortedFormula, special, data = data)

    ucTermInd <- attr (Terms, "specials")$uc  # special indices in original formula (beginning with 1 = response!)
    nUcGroups <- length (ucTermInd)
    bfpTermInd <- attr (Terms, "specials")$bfp
    nFps <- length (bfpTermInd)
    fixTermInd <- seq_along(tempVarNames)[-c(1, ucTermInd, bfpTermInd)] # indices of fixed terms
    nFixGroups <- length(fixTermInd)
    
    ## check if bfp's are present
    if (nFps == 0)
        warning(simpleWarning("no fractional polynomial terms in formula"))
    
    ## get vector with covariate entries
    vars <- attr (Terms, "variables")    # language object
    varlist <- eval (vars, envir = data)              # list
    covariates <- paste(as.list (vars)[-c(1,2)]) # vector with covariate entries (no list or response or Intercept)

    ## remove bfp() from entries and save the inner arguments
    bfpInner <- varlist[bfpTermInd]     # saved for later use
    covariates[bfpTermInd - 1] <- unlist(bfpInner) # remove bfp( ) from formula; -1 because of reponse column

    ## if ucs are present:
    if (nUcGroups){
        ## remove uc() from entries and save the inner arguments
        ucInner <- unlist(varlist[ucTermInd])
        covariates[ucTermInd - 1] <- ucInner
        
        ## determine association of terms with uc groups
        ucTermLengths <- sapply (ucInner, function (oneUc)
                                 length (attr (terms (as.formula (paste ("~", oneUc))), "term.labels"))
                                 )
        ucTermLengthsCum <- c(0, cumsum (ucTermLengths - 1)) # how much longer than 1, accumulated
        ucTermList <- lapply (seq (along = ucTermInd), function (i) # list for association uc group and assign index
                              as.integer(
                                         ucTermInd[i] - 1 + # Starting assign index
                                         ucTermLengthsCum[i]+ # add lengths from before
                                         0:(ucTermLengths[i]-1) # range for this uc term
                                         )
                              )
    } else {
        ucInner <- ucTermList <- NULL
    }
    ## consistency check:
    stopifnot(identical(length(ucTermList), nUcGroups))

    # if fixed variables are present
    if (nFixGroups) {
        ## remove uc() from entries and save the inner arguments
        fixInner <- fixTempVarNames
        covariates[fixTermInd - 1] <- fixInner

        ## determine association of terms with fixed groups
        fixTermLengths <- sapply(fixInner, function(oneFix)
            length(attr(terms(as.formula(paste("~", oneFix))), "term.labels"))
                                 )
        fixTermLengthsCum <- c(0, cumsum(fixTermLengths - 1)) # how much longer than 1, accumulated
        fixTermList <- lapply(seq(along = fixTermInd), function(i) # list for association uc group and assign index
                              as.integer(
                                         fixTermInd[i] - 1 + # Starting assign index
                                         fixTermLengthsCum[i] + # add lengths from before
                                         0:(fixTermLengths[i] - 1) # range for this uc term
                                         )
                              )
    } else {
        fixInner <- fixTermList <- NULL
    }


    ## check that not two entries are present for any covariate,
    ## i.e. the covariate names must be unique
    if(! identical(covariates,
                   unique(covariates)))
    {
        stop(simpleError(paste("Duplicate covariates in formula:",
                               paste(covariates[duplicated(covariates)],
                                     sep=", "))))
    }
    
    ## build new formula from the cleaned covariate entries
    newFormula <-                       # is saved for predict method at the end
        update (sortedFormula,
                paste (".~ 1 +",  
                       paste (covariates, collapse = "+")) # only update RHS
                )
    newTerms <- if (missing(data))
        terms(newFormula)
    else
        terms(newFormula, data = data)

    ## build model frame
    m$formula <- newTerms
    m$scale <- m$family <- m$verbose <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    ## build design matrix
    X <- model.matrix (newTerms, m)
    Xcentered <- scale(X, center=centerX, scale=FALSE)

    ## get and check weights
    weights <- as.vector(model.weights(m))
    if(is.null(weights))
        weights <- rep(1, nrow(X))
    
    if (!is.null(weights) && !is.numeric(weights)) 
        stop(simpleError("'weights' must be a numeric vector"))
    if (!is.null(weights) && any(weights < 0)) 
        stop(simpleError("negative weights not allowed"))
    
    ## get response
    Y <- model.response(m)

    ## get offsets
    offset <- as.vector(model.offset(m))
    if (!is.null(offset))
    {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    } else {
        offset <- rep.int(0, nrow(X))
    }
    
    ## check length of censoring vector in the Cox case
    if(! doGlm)
    {   
      if(class(attr(m,"na.action"))=="omit"){
        censInd <- censInd[-attr(m,"na.action")]
      }
        stopifnot(identical(length(censInd), length(Y)))
    }

    ## initialize (here e.g. the binomial matrix Y case is handled as in 'glm')
    init <- family$init(y=Y, weights=weights)

    Y <- init$y
    family$weights <- as.double(init$weights)
    family$offsets <- as.double(offset)
    family$dispersions <- as.double(family$phi / init$weights) # and zero weights ?!
    family$linPredStart <- as.double(init$linPredStart)
    
    ## which terms gave rise to which columns?
    ## (0 = intercept, 1 = first term)
    termNumbers <- attr (X, "assign") 

    ## vector of length col (X) giving uc group indices or 0 (no uc)
    ## for associating uc groups with model matrix columns: ucIndices
    fixIndices <- ucIndices <- fpMaxs <- integer (length (termNumbers))

    ## list for mapping group -> columns in model matrix: ucColList
    if(nUcGroups){
        for (i in seq (along=ucTermList)){
            ucIndices[termNumbers %in% ucTermList[[i]]] <-  i
        }

        ucColList <- lapply (seq (along=ucTermList), function (ucGroup) which (ucIndices == ucGroup))
    } else {
        ucColList <- NULL
    }

    ## list for mapping group -> columns in model matrix: fixColList
    if (nFixGroups) {
        for (i in seq(along = fixTermList)) {
            fixIndices[termNumbers %in% fixTermList[[i]]] <- i
        }

        fixColList <- lapply(seq(along = fixTermList), function(fixGroup) which(fixIndices == fixGroup))
    } else {
        fixColList <- NULL
    }


    ## vectors of length col (X) giving maximum fp degrees or 0 (no bfp)
    ## and if scaling is wanted (1) or not (0)
    ## for associating bfps with model matrix columns
    ## In addition, scale Columns or exit if non-positive values occur
    bfpInds <- bfpTermInd - 1 + attr (Terms, "intercept")           # now indexing matrix column
    for (i in seq (along=bfpInner)){
        colInd <- bfpInds[i]
        fpObj <- bfpInner[[i]]

        fpMaxs[colInd] <- attr (fpObj, "max")

        ## get scaling info
        scaleResult <- fpScale (c(attr(fpObj, "rangeVals"), # extra values not in the data
                                  X[,colInd]),              # covariate data
                                scaling = attr(fpObj, "scale")) # scaling wished?
        attr (bfpInner[[i]], "prescalingValues") <- scaleResult
        
        ## do the scaling
        X[,colInd] <- X[,colInd] + scaleResult$shift
        X[,colInd] <- X[,colInd] / scaleResult$scale

        ## check positivity
        if (min(X[,colInd]) <= 0) 
            stop (simpleError(paste("prescaling necessary for negative values in variable", fpObj)))
    }

    ## check that all maximum FP degrees are equal, so that the SWITCH move
    ## in the model sampling algorithm will always be possible.
    ## This assumption could potentially be removed later on, or otherwise the "max" option in bfp()
    ## could be removed.
    if (length(unique(fpMaxs[fpMaxs != 0])) > 1L)
        stop(simpleError("all maximum FP degrees must be identical"))

    ## Check if factor variables are used in bfp
    if(any(c('logical','factor')%in% lapply(m, class)[bfpInds]))
      stop (simpleError("logical or factor variables cannot be fractional polynomials"))
    
    
    ## check that only the intercept (one column) is a fixed term
    if (sum(! ((fpMaxs | ucIndices)| fixIndices)) > 1)
        stop (simpleError("only the intercept can be a fixed term"))
    
    ## attach a loglik-function, which is then called from the C++ code.
    ## It gives then the loglikelihood of the mu vector of means.
    family$loglik <- function(mu)
    {
        return(- 0.5 * sum(family$dev.resids(y=Y, mu=mu, wt=weights)) / phi)
    }
    ## note that this does not include normalizing constants of
    ## the sampling density, e.g. - 0.5 * log(2 * pi * phi) is *not*
    ## included in the Gaussian case!

    ## fit the null model
    nullModelFit <- glm(Y ~ 1,
                        weights=family$weights,
                        family=family$family,
                        offset=family$offsets)
    nullModelLogMargLik <-
        if(tbf)
        {
            ## if we use TBF, then the log marginal likelihood is the log
            ## test-based Bayes factor of the model versus the null model. Of
            ## course the Bayes factor of the null model versus itself is 1, so
            ## the log BF is 0.
            0
        } else {
            ## compute log marginal likelihood in the null model
            ## via the Laplace approximation:
            - nullModelFit$deviance / 2 + 0.5 * log(2 * pi * vcov(nullModelFit))  
        }
    nullModelDeviance <- nullModelFit$deviance
    
    ## compute and print cardinality of the model space to guide decision
    fpSetCards <- ifelse (fpMaxs[fpMaxs != 0] <= 3, 8, 5 + fpMaxs[fpMaxs != 0])
    getNumberPossibleFps <- function (  # computes number of possible univariate fps (including omission)
                                      maxDegree # maximum fp degree
                                      ){
        s <- ifelse (maxDegree <= 3, 8, 5 + maxDegree) # card. of power set
        singleDegreeNumbers <- sapply (0:maxDegree, function (m)
                                       choose (s - 1 + m, m))
        return (sum (singleDegreeNumbers))
    }
    singleNumbers <- sapply (fpMaxs, getNumberPossibleFps)
    totalNumber <- prod (singleNumbers) * 2^(nUcGroups) + 1 # maximum number of possible models (+1 for model with only fixed covariates)

    
    ## process the nModels argument
    if(missing(nModels))
    {
        ## then we would like to have the default number of models:
        nModels <- max(1L, floor(totalNumber / 100))
    }
    ## check nModels is at least 1
    if(nModels < 1)
    {
        stop(simpleError("nModels must at least be 1"))
    }
        
        
    ## decide if we are going to do model sampling or an exhaustive search
    if (identical(method, "ask")){
        cat ("The cardinality of the model space is at most ", totalNumber, ".\n", sep = "")
        decision <- substr (readline(paste("Do you want to do a deterministic search for the best model (y)",
                                           "or sample from the model space (n) or abort (else) ?\n")),
                            1, 1)

        ## ensure that correct decision string has been entered, otherwise abort
        if(! decision %in% c("y", "n"))
        {
            cat("Aborting.\n")
            return()
        }

    } else {
        decision <- switch(method, exhaustive = "y", sampling = "n")
    }
    ## and the match.arg before ensures that no further problems can occur.

    ## translate the decision to logical variable 
    doSampling <- identical(decision, "n")
    
    ## ensure that we only do model sampling if there is more than 1 FP term in the model
    if(doSampling && identical(nFps, 1L))
    {
        warning(simpleWarning(paste("We need to do an exhaustive computation of all models,",
                                    "because there is only 1 FP term in the model!")))
        doSampling <- FALSE
    }

    ## if sampling, we possibly ask for the chainlength
    if (doSampling)
    {
        ## get chainlength?
        if (identical(method, "ask"))
        {
            chainlength <- as.numeric (readline ("How long do you want the Markov chain to run?\n"))
        }
      
        
        ## compute the default number of models to be saved
        if(is.null(nModels))
        {
            nModels <- as.integer(max(chainlength / 100, 1L))
        }
        else
            stopifnot(nModels >= 1L)

        ## check the chosen cache size
        nCache <- as.integer(nCache)
        stopifnot(nCache >= nModels)        
        
        if (verbose)
        {
            cat("Starting sampler...\n")            
        }
    }
    else                                
    {
        if (verbose)
        {
            cat("Starting with computation of every model...\n")
        }
    }

    ## start the progress bar (is continued in the C++ code)
    if(verbose)
    {
        cat ("0%", rep ("_", 100 - 6), "100%\n", sep = "")
    }

    ## if Cox model is requested,
    ## order the data according to the survival times!
    if(! doGlm)
    {
        sorted <- order(Y)
        
        Y <- Y[sorted]
        censInd <- censInd[sorted]
        X <- X[sorted, ]
        Xcentered <- Xcentered[sorted, ]

        if(! all(sorted == seq_along(Y)))
        {
            warning("Input data were reordered so that the survival times are sorted") 
        }
    }

    ## pack the data together
    data <- list(x=X,                   # design matrix
                 xCentered=Xcentered,   # centered design matrix
                 y=as.double(Y),        # response vector
                 nObs=nrow(X),          # number of observations
                 censInd=as.integer(censInd)) # binary censoring indicator

    ## pack the FP info things together
    fpInfos <- list(fpmaxs=as.integer(fpMaxs[fpMaxs != 0]), # vector of maximum fp degrees)
                    fppos=as.integer(bfpInds),              # vector of fp columns
                    fpcards=as.integer(fpSetCards), # cardinality of corresponding power sets
                    fpnames=if (nFps == 0) character(0) else unlist(bfpInner)) # names of fp terms. Note that
                                        # it is necessary to have this if-else
                                        # construction (ifelse does not work,
                                        # because it would take the first
                                        # element of the length-0 vector
                                        # character(0), which is NA!!)
    
    
    
    ## pack the UC info things together
    ucInfos <- list(ucIndices=as.integer(ucIndices), # vector giving uncertainty custer indices
                                        # (column -> which group) 
                    ucColList=ucColList) # list for group -> which columns mapping


    ## pack the UC info things together
    fixInfos <- list(fixIndices = as.integer(fixIndices), # vector giving fixed covaraite custer indices
                                    # (column -> which group) 
                fixColList = fixColList) # list for group -> which columns mapping



    ## pack model search configuration:
    searchConfig <- list(totalNumber=as.double(totalNumber), # cardinality of model space                         
                         nModels=as.integer(nModels),         # number of best
                                        # models returned
                         empiricalBayes=empiricalBayes, # use EB for g and
                                        # conditional marginal likelihoods?
                         useFixedg=useFixedg,   # use a fixed value of g?
                         useFixedc = fixedcfactor,
                         doSampling=doSampling,               # shall model sampling be done? If
                                        # false, then exhaustive search.
                         chainlength=as.double(chainlength),  # how many times should a jump be
                                        # proposed?
                         nCache=nCache, # how many models to cache at the same time
                         largeVariance=as.double(largeVariance), # what is a "large" variance output
                                        # of BFGS?
                         useBfgs=useBfgs) # should we use the BFGS algorithm (or
                                        # Brent's optimize)?                        

    ## pack prior and likelihood information:
    distribution <- list(nullModelLogMargLik=nullModelLogMargLik, # log marg lik
                                        # of the null model.
                         nullModelDeviance=nullModelDeviance, # corresponding deviance
                         doGlm=doGlm,   # should we do a GLM or a Cox model search?
                         tbf=tbf,       # should TBF methodology be used to
                                        # compute Bayes factors?
                         fixedg=as.double(fixedg), # the fixed value of g
                                        # (0 if not used)
                         gPrior=priorSpecs$gPrior, # prior on the covariance
                                        # factor g (S4 class object)
                         modelPrior=priorSpecs$modelPrior, # model prior string                         
                         family=family,    # GLM family and link,
                         yMean=mean(Y),   #  pass the mean, which we use when fixedc=TRUE
                         empiricalgPrior=empiricalgPrior) # should we use empirical g prior

    ## pack other options
    options <- list(verbose=verbose,           # should progress be displayed?
                    debug=debug,               # echo debug-style messages?
                    gaussHermite=gaussHermite,   # nodes and weights for Gauss
                                        # Hermite quadratures
                    useOpenMP=useOpenMP, # should we use openMP for speed up?
                    higherOrderCorrection=higherOrderCorrection) # should
                                        # the higher-order Laplace correction be used?    
    
    ## then go C++
    Ret <- cpp_glmBayesMfp(data,
                           fpInfos,
                           ucInfos,
                           fixInfos,
                           searchConfig,
                           distribution,
                           options)

    ## C++ attaches the following attributes:

    ## numVisited
    ## inclusionProbs
    ## logNormConst

    ## name the inclusion probabilities
    names (attr (Ret, "inclusionProbs")) <- c(unlist (bfpInner), ucInner)

    ## name the models with the model index
    names (Ret) <- 1:length(Ret)

    ## attach additional information:

    ## information passed to C++, which is important for e.g. the function "sampleGlm"
    attr (Ret, "data") <- data
    attr(Ret, "fpInfos") <- fpInfos
    attr(Ret, "ucInfos") <- ucInfos
    attr(Ret, "fixInfos") <- fixInfos
    attr(Ret, "searchConfig") <- searchConfig
    attr(Ret, "distribution") <- distribution
    attr(Ret, "options") <- options

    ## original call and formula
    attr (Ret, "call") <- call
    attr (Ret, "formula") <- newFormula

    ## prior specs argument
    attr (Ret, "priorSpecs") <- priorSpecs


        ## list with index info
    #fixedInds <- setdiff (1:ncol (X), c (bfpInds, which (ucIndices > 0)))
    attr (Ret, "indices") <- list (uc = ucIndices,
                                   ucList = ucColList,
                                   bfp = bfpInds,
                                   fixed = fixIndices)

    
    ## names of the terms
    fixedNamesInds <- c(setdiff(2:length (varlist), unlist (attr (Terms, "specials"))) , 1)

    # interceptName <- ifelse (attr (Terms, "intercept"), "(Intercept)", NULL)
    if (doGlm) { interceptName <- "(Intercept)"
    } else if (!doGlm) interceptName <- NULL
    
    
    attr (Ret, "termNames") <- list (fixed=c(interceptName, fixInner),
                                     bfp=unlist(bfpInner),
                                     uc=ucInner)

    ## matrix with shift/scale info, maximum degree and cardinality of powerset
    shiftScaleMaxMat <- matrix (nrow = nFps, ncol = 4)
    colnames (shiftScaleMaxMat) <- c ("shift", "scale", "maxDegree",
                                      "cardPowerset")

    if(nFps > 0L)
    {
        shiftScaleMaxMat[, 1:2] <- matrix(unlist(lapply(bfpInner,
                                                        attr,
                                                        "prescalingValues")),
                                          ncol = 2,
                                          byrow = TRUE)
        shiftScaleMaxMat[, 3] <- fpMaxs[fpMaxs != 0]
        shiftScaleMaxMat[, 4] <- fpSetCards

        rownames (shiftScaleMaxMat) <- unlist (bfpInner)
    }

    attr (Ret, "shiftScaleMax") <- shiftScaleMaxMat

    ## set class and return
    class (Ret) <- c("GlmBayesMfp", "list")
    return(Ret)
}


