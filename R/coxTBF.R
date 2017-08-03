#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs for GLMs
## 
## Description:
## Friendlier interface to fit Cox models with glmBayesMfp
##
## History:
## 14/07/2015 Copy from CoxTBFs project
#####################################################################################

##' @include helpers.R
##' @include getDesignMatrix.R
{}


##' Fit Cox models using glmBayesMfp
##'
##' A simplified formula based interface to \code{\link{glmBayesMfp}} to fit Cox models. Can return 
##' Maximum a posteriori (MAP) model, Median probability model (MPM) or Bayesian model average (BMA).
##' Provides global empirical Bayes and AIC/BIC based model inference.
##' 
##' @param formula model formula with Surv object as LHS and \code{\link{uc}} or \code{\link{bfp}} 
##' variables as RHS.
##' @param data data.frame for model variables
##' @param type type of model to fit, one of "MAP","MPM","BMA","BMAFull"
##' @param baseline how to calculate the baseline hazard function. "cox" uses unshrunken
##' coefficients. "shrunk" refits baseline with shrunken coefficients (default).
##' @param globalEB use global empirical bayes estimate of g (default=FALSE)
##' @param IC use information criteria based model selection (default=FALSE). Either "AIC" or "BIC".
##' @param sep estimate baseline hazard for each estimate of model coefficients (default=FALSE).
##' @param keepModelList keep the model list returned by glmBayesMfp for MAP and MPM models (default=FALSE).
##' @param ... additional arguments to pass to \code{\link{glmBayesMfp}}
##' @param overrideConfig replaces the the MAP model with the given configuration, which is passed to \code{\link{computeModels}}
##'
##' @return An object of S3 class \code{TBFcox} or \code{TBFcox.sep} if sep=TRUE.
##' 
##' @keywords models regression
##' 
##' @import stats
##' @export

coxTBF <- function(formula, data, type, baseline='shrunk', globalEB=FALSE, IC=FALSE, sep=FALSE, keepModelList = FALSE, ..., overrideConfig){
 
  formula <- as.formula(formula)
  
  LHS <- formula[[2]][[2]]
  RHS <- paste(attr(terms(formula),"term.labels"),collapse=" + ")
  selection.formula  <- formula(paste(LHS,"~",RHS))
  
  tryCatch(type<-match.arg(type,c("MAP","MPM","BMA","BMAfull")), error= function(e) print("Invalid type specification"))
  tryCatch(baseline<-match.arg(baseline,c("cox","shrunk")), error= function(e) print("Invalid baseline specification"))
  
  time.var = as.character(formula[[2]][[2]])
  status.var = as.character(formula[[2]][[3]])
  
  #Set up the g-prior
  nEvents<- sum(data[[status.var]])
  
  prior.hypergn <- CustomGPrior(logDens = function(g) {
    return(-log(nEvents) - 2 * log1p(g/nEvents))
  })
  
  #############################################################################################
  #Handle Global Empirical Bayes
  #First we need to fit the model, then we can hand a fixed g to the normal search procedure
  if(globalEB==TRUE){
    
    full.model.space <- 2^length(labels(terms(selection.formula)))
    
    gEB.models <- glmBayesMfp(selection.formula,
                              censInd = data[[status.var]],
                              data = data,
                              tbf = TRUE,
                              empiricalBayes=TRUE,
                              priorSpecs = list(gPrior = prior.hypergn, modelPrior = "dependent"),
                              nModels = full.model.space,
                              method = "exhaustive",
                              verbose = FALSE)
    
    k <- length(gEB.models)
    deviances <- sapply(1:k, function(i) gEB.models[[i]]$information$residualDeviance)
    log.prior.prob <- sapply(1:k, function(i) gEB.models[[i]]$information$logPrior)
    
    ucList <- attr(gEB.models, "indices")$ucList
    degrees <- sapply(1:k, function(i)
      length(unlist(ucList[gEB.models[[i]]$configuration$ucTerms])))
    valid <- is.na(deviances)!=TRUE
    
    TBF <- function(g) sum( ((g+1)^(-degrees/2) * exp(g/(g+1) * deviances/2)*exp(log.prior.prob))[valid])
    
    globalEB.fun <- function(g) sapply(g,TBF)
    
    bestg <- optimise(f = globalEB.fun, interval=c(0.5,1000), maximum=TRUE)$maximum
    print(paste("Global EB chooses g =",bestg))
  }
  #end globalEB
  #############################################################################################
  
  args1 <- list(selection.formula,
                censInd = data[[status.var]],
                data = data,
                tbf = TRUE,
                priorSpecs = list(gPrior = prior.hypergn, modelPrior = "dependent"),
                chainlength = 5000,
                method = "sampling",
                verbose = FALSE)
  
  inargs <- list(...)
  
  if(IC !=FALSE) inargs <- c(inargs,fixedg=10*nrow(data))
  
  if(globalEB && exists("bestg")){
    inargs <- c(inargs,fixedg=bestg)
  }
  
  #match args and ...
  if(length(inargs)){
    for(i in 1:length(inargs))
      args1[[names(inargs)[i]]] <- inargs[[i]]
  }
  
  #fit the model using TBF method
  models <- do.call(glmBayesMfp, args1)
  
  this.call <- match.call()
  
  if(globalEB) {
    geb <- which(names(this.call)=="globalEB")
    names(this.call)[geb] <- "fixedg"
    this.call[geb] <- bestg
  }
  
  # If we are doing AIC/BIC we need to fit models, calculate IC values and new probabilities.
  if(IC != FALSE & type != "BMAfull"){ 
    IC.values <- numeric(length(models))
    nullModel <- NA
    
    
    design.names <- colnames(attr(models,"data")$xCentered)
    
    new.data <- data.frame(cbind(status.var=attr(models,"data")$censInd,
                                 time.var=attr(models,"data")$y,
                                 attr(models,"data")$xCentered[,-1]))
    
    colnames(new.data)[1:2] <- c(status.var,time.var)
    
    
    
    for(j in 1:length(models)){
      # new.design.matrix <- getDesignMatrix(object=models[j])[,-1,drop=FALSE]
      # 
      # new.data <- data.frame(cbind(status.var=attributes(models[j])$data$censInd,
      #                              time.var=attributes(models[j])$data$y,
      #                              new.design.matrix))
      # 
      # colnames(new.data)[1:2] <- c(status.var,time.var)
      
      #covariate indices corresponding to data frame
      ind <- unlist(attr(models,"indices")$ucList[models[[j]]$configuration$ucTerms])
      
      
      
      model.formula <- paste("survival::Surv(",time.var,",",status.var,")~1",paste(c(" ", design.names[ind]), collapse = " + "))
      
      model.cph <- rms::cph(formula(model.formula), data=new.data, surv=FALSE, se.fit = FALSE, y=FALSE, x=FALSE)
      
      # if(j%%10==0) print(j)
     
      
      if(IC=='AIC') IC.values[j] <- -AIC(model.cph)/2
      
      if(IC=='BIC') IC.values[j] <- -AIC(model.cph, k = log(sum(new.data[status.var])))/2
      
      #save this model so we can delete it
      if(is.na(IC.values[j])| length(models[[j]]$configuration$ucTerms)==0) {
        nullModel <- j
      }
      
    }
    
    #delete null model since it is can't be used later
    if(!is.na(nullModel)){
      IC.values <- IC.values[-nullModel]
      models <- models[-nullModel]
    }
    models.df <- as.data.frame(models)
    
    
    IC.posteriors <- exp((IC.values-max(IC.values))+models.df$logPrior)
    IC.posteriors <- IC.posteriors/sum(IC.posteriors)
    
    if(type=="MPM"){
      attr(models, "inclusionProbs") <- apply(models.df[,-(1:3)], 2, function(x) sum(x*IC.posteriors))
    }
    
    if(type=="MAP"){
      models[1] <- models[which.max(IC.posteriors)]
      #also include inclusion probs for *IC methods
      attr(models, "inclusionProbs") <- apply(models.df[,-(1:3)], 2, function(x) sum(x*IC.posteriors))
    }
  }
  #End AIC/BIC
  

  
  
    
  if(type=="MAP"){
    model.listpart <- models[1]

    #if overrideConfig exists, we use it instead of MPM to fit the desired model.
    if(!missing(overrideConfig)) {
      print("Using overrideConfig model")
      
      new.models <- computeModels(list(overrideConfig), models)
      model.listpart <- new.models[1]
    }
    
  } else if(type=="MPM"){
    #what are the covariates with inclusion probability >0.5?
    mpm.vars <- attr(models, "inclusionProbs")>0.5
    #which model includes only these?
    model.df <- as.data.frame(models)[,-c(1:3)]
    mpm.index <-which(sapply(1:length(models), function(i) all(model.df[i,] == mpm.vars)))
    if(!any(mpm.index)) {
      #if the mpm wasn't fitted we must construct it
      print("MPM model wasn't fitted so we construct it.")
      mpm.model.config <- models[[1]]$configuration
      mpm.model.config$ucTerms <- which(mpm.vars)

      new.models <- computeModels(list(mpm.model.config), models)
      model.listpart <- new.models[1]
      print("Success.s")
    } else {
      model.listpart <- models[mpm.index]
    }
    
  } else if(type=="BMAfull"){
    # BMA-FULL ----------------------------------------------------------------   
    bma.surv <- list()
    bma.coefs <- list()
    bma.survfits <- list()
    surv.time <- numeric()
    
    
    
    models[ which(sapply(1:length(models),
                         function(i) length(models[i][[1]]$configuration$ucTerms))==0)] <- NULL
    
    #if we are doing AIC or BIC we need to calculate (and store) the values separately
    if(IC!=FALSE) IC.values <- numeric(length(models))
    
    for(j in 1:length(models)){
      #    foreach(j = 1:length(models)) %dopar% {
      new.design.matrix <- getDesignMatrix(object=models[j], intercept=FALSE)
      
      new.data <- data.frame(cbind(status.var=attributes(models[j])$data$censInd,
                                   time.var=attributes(models[j])$data$y,
                                   new.design.matrix))
      
      colnames(new.data)[1:2] <- c(status.var,time.var)
      model.formula <- paste("survival::Surv(",time.var,",",status.var,")~.")
      
      #calculate values for each model
      bma.coefs[[j]] <- getModelCoefs(models[j], 
                                        mcmc=McmcOptions(burnin=0L, step=1L, samples=100))
      # print(colnames(new.data))
      model.cph <- rms::cph(formula(model.formula), data=new.data, surv=TRUE, se.fit = FALSE, y=TRUE, x=TRUE)
      if(j==1) surv.time <- model.cph$time
      
      if(baseline=="cox") {
        bma.surv[[j]] <- rms::Survival(model.cph)
        
        shrunk.survfit <- survival::survfit(model.cph)
        shrunk.surv <- c(1,shrunk.survfit$surv)
        if(length(shrunk.surv) < length(model.cph$surv)) shrunk.surv<- c(shrunk.surv, shrunk.surv[length(shrunk.surv)])
        bma.survfits[[j]] <- shrunk.surv
      }
      
      if(baseline=="shrunk"){
        
        # Take cox model object, put in shrunken coefficients, re-estimate baseline hazard
        shrunk.cph <- model.cph
        shrunk.mm <- new.design.matrix  #model.matrix(model.formula, data=data)[,-1] 
        #     shrunk.mm <-scale(shrunk.mm)
        shrunk.cph$linear.predictors <- shrunk.mm %*% bma.coefs[[j]]
        
        
        # shrunk.survfit <- survival::survfit(shrunk.cph)
        
        shrunk.survfit <- try(survival::survfit(shrunk.cph))
        
        if(class(shrunk.survfit)[1]=="try-error"){
          this.model <- models[j]
          these.coefs <- bma.coefs[[j]]
          save("model.cph","shrunk.mm","these.coefs","this.model",file="Faulty.Files.Rdata")
          
          g <- exp(models[[j]]$information$zMode)
          shrunk.cph$linear.predictors <- shrunk.mm %*% (model.cph$coefficients * g/(g+1) )
          shrunk.survfit <- try(survival::survfit(shrunk.cph))
          
        }        
        
        
        shrunk.surv <- c(1,shrunk.survfit$surv)
        if(length(shrunk.surv) < length(model.cph$surv)) shrunk.surv<- c(shrunk.surv, shrunk.surv[length(shrunk.surv)])
        shrunk.cph$surv <- shrunk.surv
        
        bma.surv[[j]] <- rms::Survival(shrunk.cph)
        
        bma.survfits[[j]] <- shrunk.surv
      }
      
      
      if(IC!=FALSE){
        if(IC=='AIC') IC.values[j] <- -AIC(model.cph)/2
        if(IC=='BIC') IC.values[j] <- -AIC(model.cph, k = log(sum(new.data[status.var])))/2
        
      }
      
    }
    
    #unlist gives the factor levels the wrong names, so take the better ones
    #from cph to avoid confusion!
    #   names(model.coefs) <- names(model.cph$coefs)
    
    ret <- list()
    
    ret$formula <- writeFormula(models[1], time.var, status.var) #model.formula
    ret$coefs <- bma.coefs
    ret$data <- data
    ret$call <- this.call
    
    
    if(IC==FALSE) {
      ret$probability <- posteriors(models)
    } else if(IC!=FALSE){
      log.priors <- unlist(lapply(models, function(i) i$information$logPrior))
      IC.posteriors <- exp((IC.values-max(IC.values)) + log.priors)
      ret$probability <- IC.posteriors/sum(IC.posteriors)
    }
    
    ret$survfit <- matrix(unlist(bma.survfits), byrow=TRUE, nrow=length(bma.survfits))
    ret$time <- surv.time
    ret$bma.surv <- bma.surv
    ret$model.object <- models #just save the first one for space
    class(ret) <- "TBFcox.BMA"
    return(ret)
    
    
  }  else if(type=="BMA"){
    
    
    #check for null models and remove them. The fitting procedure can't handle them
    #this might be a bit dodgy, but no alternative solution comes to mind
    isNull <- apply(as.data.frame(models)[,-(1:3)],1,sum)==0
    models[isNull] <- NULL
    
    if(IC != FALSE) {
      sbma <- sampleBma(models, postProbs=IC.posteriors[!isNull])
    } else {
      sbma <- sampleBma(models)
    }
    print("1")
    
    #     uc.included <- which(unlist(lapply(sbma$samples@ucCoefs, function(x) is.numeric(mean(x)))))
    uc.included <- which(apply(as.data.frame(models)[,-c(1:3)],2,any))
    
    
    fake.model <- models[1]
    fake.model[[1]]$information <- list()
    
    uc.num.inc <- which(attr(fake.model,"termNames")$uc %in% names(uc.included))
    
    fake.model[[1]]$configuration <- list(powers=list(), ucTerms=uc.num.inc)
    
    new.design.matrix <- getDesignMatrix(object=fake.model, intercept=FALSE)
    
    
    #get coefficients for BMA
    #     if(BMA.method=="coefs"){
    #       bma.coefs <- unlist(lapply(sbma$samples@ucCoefs, rowSums))/sbma$samples@nSamples
    #     }
    
    #     if(BMA.method=="fits"){
    bma.coefs <- .lm.fit(x=new.design.matrix, y=rowMeans(sbma$samples@fitted))$coefficients
    #     }
    
    #things were out of order making predictions terrible!
    #     bma.coefs <- bma.coefs[order(names(bma.coefs))]
    names(bma.coefs) <- colnames(new.design.matrix)
    
    new.data <- data.frame(cbind(status.var=attributes(fake.model)$data$censInd,
                                 time.var=attributes(fake.model)$data$y,
                                 new.design.matrix))
    
    colnames(new.data)[1:2] <- c(status.var,time.var)
    model.formula <- paste("survival::Surv(",time.var,",",status.var,")~.")
    
    model.cph <- rms::cph(formula(model.formula), data=new.data, surv=TRUE, se.fit = FALSE, y=TRUE, x=TRUE)
    
    if(baseline=="cox") {
      bma.surv <- rms::Survival(model.cph)
    }
    
    if(baseline=="shrunk"){
      # Take cox model object, put in shrunken coefficients, re-estimate baseline hazard
      shrunk.cph <- model.cph
      shrunk.mm <- new.design.matrix  #model.matrix(model.formula, data=data)[,-1] 
      #     shrunk.mm <-scale(shrunk.mm)
      shrunk.cph$linear.predictors <- shrunk.mm %*% bma.coefs
      shrunk.survfit <- survival::survfit(shrunk.cph)
      shrunk.surv <- c(1,shrunk.survfit$surv)
      if(length(shrunk.surv) < length(model.cph$surv)) shrunk.surv<- c(shrunk.surv, shrunk.surv[length(shrunk.surv)])
      shrunk.cph$surv <- shrunk.surv
      
      bma.surv <- rms::Survival(shrunk.cph)
    }
    
    #unlist gives the factor levels the wrong names, so take the better ones
    #from cph to avoid confusion!
    #   names(model.coefs) <- names(model.cph$coefs)
    ret <- list()
    
    #ret$formula <- writeFormula(models[1], time.var, status.var) #model.formula
    ret$formula <- formula(paste("survival::Surv(",time.var,",", status.var,") ~", paste(paste(names(sbma$samples@ucCoefs)),collapse=" + ")))
    ret$coefs <- bma.coefs
    ret$data <- data
    ret$call <- this.call
    
    
    ret$survival <- bma.surv
    ret$model.object <- fake.model #just save the first one for space
    class(ret) <- "TBFcox"
    return(ret)
  }
  #############################################################################################################
  
  if(type %in% c("MAP","MPM")){
    if(!keepModelList) rm(models)
    
    new.design.matrix <- getDesignMatrix(object=model.listpart, intercept=FALSE)
  
    new.data <- data.frame(cbind(status.var=attributes(model.listpart)$data$censInd,
                                 time.var=attributes(model.listpart)$data$y,
                                 new.design.matrix))
    
    colnames(new.data)[1:2] <- c(status.var,time.var)
    model.formula <- paste("survival::Surv(",time.var,",",status.var,")~.")
    
    #Do the original shortcut way, using E(beta) in the exp.
    if(sep==FALSE){
       model.coefs <- getModelCoefs(model.listpart, sep=FALSE)
       model.cph <- rms::cph(formula(model.formula), data=new.data, surv=TRUE, se.fit = FALSE, y=TRUE, x=TRUE)
    
       ret <- list()
    
       ret$formula <- writeFormula(model.listpart, time.var, status.var) #model.formula
       ret$coefs <- model.coefs
       ret$data <- data
       ret$call <- this.call
       ret$survival <- function(){}
       ret$model.object <- model.listpart
       if(exists("models")) ret$model.list <- models
    
      if(baseline=="shrunk"){
        # Take cox model object, put in shrunken coefficients, re-estimate baseline hazard
        shrunk.cph <- model.cph
        shrunk.mm <- new.design.matrix  #model.matrix(model.formula, data=data)[,-1] 

        shrunk.cph$linear.predictors <- shrunk.mm %*% model.coefs
        shrunk.survfit <- survival::survfit(shrunk.cph)
        shrunk.surv <- c(1,shrunk.survfit$surv)
        if(length(shrunk.surv) < length(model.cph$surv)) shrunk.surv<- c(shrunk.surv, shrunk.surv[length(shrunk.surv)])
        shrunk.cph$surv <- shrunk.surv
      
        ret$survival <- rms::Survival(shrunk.cph)
      
      } else if(baseline=="cox"){
        ret$survival <- rms::Survival(model.cph)
      }
      
      class(ret) <- "TBFcox"
      return(ret)
    }
    
    #Do the full way with E(h_i(t)exp(beta))
    if(sep==TRUE){
      model.coefs <- getModelCoefs(model.listpart, sep=TRUE)
      model.cph <- rms::cph(formula(model.formula), data=new.data, surv=TRUE, se.fit = FALSE, y=TRUE, x=TRUE)
            
      ret <- list()
      
      ret$formula <- writeFormula(model.listpart, time.var, status.var) #model.formula
     
      ret$data <- data
      ret$call <- this.call
      
      ret$model.object <- model.listpart
      if(exists("models")) ret$model.list <- models
      
      ret$survival <- vector(mode="list", length=length(model.coefs))
      ret$coefs <- model.coefs
       
              
        if(baseline=="shrunk"){
          # Take cox model object, put in shrunken coefficients, re-estimate baseline hazard
          shrunk.cph <- model.cph
          shrunk.mm <- new.design.matrix  #model.matrix(model.formula, data=data)[,-1] 
          
          for(k in 1:length(model.coefs)){
            shrunk.cph$linear.predictors <- shrunk.mm %*% model.coefs[[k]]
            shrunk.survfit <- survival::survfit(shrunk.cph)
            shrunk.surv <- c(1,shrunk.survfit$surv)
            if(length(shrunk.surv) < length(model.cph$surv)) shrunk.surv <- c(shrunk.surv, shrunk.surv[length(shrunk.surv)])
            shrunk.cph$surv <- shrunk.surv
            
            ret$survival[[k]] <- rms::Survival(shrunk.cph)
          }
          
        } else if(baseline=="cox"){
          for(k in 1:length(model.coefs)){
            ret$survival[[k]] <- rms::Survival(model.cph)
          }
      }
      
      class(ret) <- "TBFcox.sep"
      return(ret)
    }
    
  
  }#end MAP/MPM
}#end function