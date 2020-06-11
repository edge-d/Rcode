################################################################################################


################################################################################################


qErr <- function(predictor, predictand, quant){
  
  #if input as data.frame - convert to vectors
  if(class(predictor) == "data.frame"){
    predictor <- predictor[,1]
  }
  if(class(predictand) == "data.frame"){
    predictand <- predictand[,1]
  }
  
  
  #Initiate variables
  sumsqErr <- 0
  
  errorList <- rep(NA, length(predictor))
  for (i in 1:length(predictor)){
    sqErr <- (predictor[i] - predictand[i])^2
    sumsqErr <- sumsqErr + sqErr
    errorList[i] <- sqrt(sqErr)
  }
  qErr <- quantile(errorList, c(quant))

  return(qErr)
}


################################################################################################


################################################################################################