#' Prediction methods for CoxTBF objects
#' 
#' Predicts survival probabilities at given times. Compatible with predictSurvProb functions 
#' required by \code{pec} package.
#'
#' @param object a model fitted with \code{\link{coxTBF}}
#' @param newdata a dataframe with the same variables as the original data used to fit the object
#' @param times a vector of times to predict survival probability for
#' @param ... not used.
#'
#' @return A data frame of survival probabilities with rows for each row of newdata and columns for each time.
#' @export
#'
#' 
predict.TBFcox <- function(object, newdata, times, ...){
  # print("predictSurvProb TBF.Cox")
  # print(paste("Train N =",nrow(attr(object$model.object,"data")$x),
  #             ". Predict N = ", nrow(newdata)))
  # From sampleGlm.R from glmBfp packages by Sabanes Bove
  obj <- object$model.object
  config <- obj[[1]]$configuration
  ## get new data matrix
  tempX <- constructNewdataMatrix(object=obj,
                                  newdata=newdata)
  
  ## copy old GlmBayesMfp object
  tempObj <- obj
  
  ## correct model matrix in tempMod to new data matrix
  attr(tempObj, "data") <-
    list(x=tempX,
         xCentered=scale(tempX, center=TRUE, scale=FALSE)) 
  
  ## so we can get the design matrix
  
  #lets see what happens if we do this another way
  newDesignUC <- getDesignMatrix(modelConfig=config,
                                           object=tempObj,
                                           intercept=FALSE,
                                           center = FALSE)
  
  
  oldDesignUC <- getDesignMatrix(modelConfig=config,
                                 object=obj,
                                 intercept=FALSE,
                                 center=FALSE)
  
  oldMeans <- colMeans(oldDesignUC)
  
  for(i in 1:length(oldMeans)){
    newDesignUC[,i] <- newDesignUC[,i] - oldMeans[i]
  }
  
  ## so the linear predictor samples are
  linPredSamples <- newDesignUC %*% object$coefs
  
  
  return(object$survival(times, linPredSamples))
}





#' Prediction methods for CoxTBF objects with separate estimates
#' 
#' Predicts survival probabilities at given times. Compatible with predictSurvProb functions 
#' required by \code{pec} package. Predicts objects with fitted with \code{sep=TRUE}
#'
#' @param object a model fitted with \code{\link{coxTBF}}
#' @param newdata a dataframe with the same variables as the original data used to fit the object
#' @param times a vector of times to predict survival probability for
#' @param ... not used.
#'
#' @return A data frame of survival probabilities with rows for each row of newdata and columns for each time.
#' @export
#'
#' 
predict.TBFcox.sep <- function(object, newdata, times, ...){
  print("predictSurvProb TBF.Cox.sep")
  print(paste("Train N =",nrow(attr(object$model.object,"data")$x),
              ". Predict N = ", nrow(newdata)))
  # From sampleGlm.R from glmBfp packages by Sabanes Bove
  obj <- object$model.object
  config <- obj[[1]]$configuration
  ## get new data matrix
  tempX <- constructNewdataMatrix(object=obj,
                                  newdata=newdata)
  
  ## copy old GlmBayesMfp object
  tempObj <- obj
  
  ## correct model matrix in tempMod to new data matrix
  attr(tempObj, "data") <-
    list(x=tempX,
         xCentered=scale(tempX, center=TRUE, scale=FALSE)) 
  
  ## so we can get the design matrix
  
  #lets see what happens if we do this another way
  newDesignUC <- getDesignMatrix(modelConfig=config,
                                           object=tempObj,
                                           intercept=FALSE,
                                           center=FALSE)
  
  
  oldDesignUC <- getDesignMatrix(modelConfig=config,
                                 object=obj,
                                 intercept=FALSE,
                                 center=FALSE)
  
  oldMeans <- colMeans(oldDesignUC)
  
  for(i in 1:length(oldMeans)){
    newDesignUC[,i] <- newDesignUC[,i] - oldMeans[i]
  }
  
  ## so the linear predictor samples are
  linPredSamples <- lapply(object$coefs, function(coefs) newDesignUC %*% coefs)
  
  k <- length(linPredSamples)
  
  ret.survival <- lapply(1:k, function(i) object$survival[[i]](times, linPredSamples[[i]]) )
  
  ret.surv2 <- matrix(0, nrow = nrow(as.data.frame(ret.survival[[1]])), ncol =  length(times))
  for(i in 1:k){
    ret.surv2 <- ret.surv2 + ret.survival[[i]]/k
  }
  
  return(ret.surv2)
}




#' Prediction methods for CoxTBF objects for BMA models
#' 
#' Predicts survival probabilities at given times. Compatible with predictSurvProb functions 
#' required by \code{pec} package. Predicts BMA objects.
#'
#' @param object a model fitted with \code{\link{coxTBF}}
#' @param newdata a dataframe with the same variables as the original data used to fit the object
#' @param times a vector of times to predict survival probability for
#' @param ... not used.
#'
#' @return A data frame of survival probabilities with rows for each row of newdata and columns for each time.
#' @export
#'
#' 
predict.TBFcox.BMA <- function(object, newdata, times, ...){
  post <- object$probability / sum(object$probability)
  time <- object$time
  
  N <- nrow(newdata)
  k <- length(object$model.object)
  t <- length(times)
  
  preds <- matrix(0, nrow=t, ncol=N)
  Big.Matrix <- matrix(0, nrow = k, ncol = t)
  LP.Matrix <- matrix(0, nrow = N, ncol = k)
  W.Matrix <- diag(x=post)
  
  
  print(paste("Predict BMA","N=",N,"k=",k,"t=",t))
  #Populate (exp) LP.Matrix 
  #This can probably be made faster since I think some of the matrices are identical
  
  print("Making Model Matrices")
  for(i in 1:length(object$model.object)){
    # From sampleGlm.R from glmBfp packages by Sabanes Bove
    obj <- object$model.object[i]
    config <- obj[[1]]$configuration
    ## get new data matrix
    tempX <- constructNewdataMatrix(object=obj,
                                    newdata=newdata)
    
    ## copy old GlmBayesMfp object
    tempObj <- obj
    
    ## correct model matrix in tempMod to new data matrix
    attr(tempObj, "data") <-
      list(x=tempX,
           xCentered=scale(tempX, center=TRUE, scale=FALSE)) 
    
    ## so we can get the design matrix
    
    #lets see what happens if we do this another way
    newDesignUC <- getDesignMatrix(modelConfig=config,
                                   object=tempObj,
                                   intercept=FALSE,
                                   center=FALSE)
    
    
    oldDesignUC <- getDesignMatrix(modelConfig=config,
                                   object=obj,
                                   intercept=FALSE,
                                   center=FALSE)
    
    oldMeans <- colMeans(oldDesignUC)
    
    newDesignUC <- newDesignUC - rep(oldMeans, each=nrow(newDesignUC))
    
    ## so the linear predictor samples (after exp transform) are
    LP.Matrix[,i] <- exp(newDesignUC %*% object$coefs[[i]])
  }
  
  print("Making Big Matrix")
  #Populate Big.Matrix
  for (i in 1:t) {
    tm <- max((1:length(time))[time <= times[i] + 1e-06])
    Big.Matrix[,i] <- object$survfit[,tm]
    if (times[i] > max(time) + 1e-06) 
      Big.Matrix[,i] <- NA
  }
  
  print("Making Predictions")
  #   for(i in 1:k){
  #     print(paste("p=",i,"of",k))
  #     preds <- preds + post[i]* outer(Big.Matrix[i,], LP.Matrix[,i], "^")
  #   }
  preds <- predBMAcpp(Big.Matrix, LP.Matrix, post)
  
  print("Finished Predictions")
  preds <- t(preds)
  if (is.data.frame(preds)) colnames(preds) <- times
  return(preds)
}

