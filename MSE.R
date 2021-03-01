
##Function to calculate MSE, RE, and CE
#Significance of RE and CE stats based on Ebisuzaki MC from:
#https://doi.org/10.1016/j.dendro.2011.08.003

#recon: reconstruction covering full period of interest
#target: reconstruction target, also covering full period of interest
#calibInt: start and end years of calibration interval, such as c(1,52)
#MC: Number of Monte Carlo iterations to calculate significance levels
#alpha: significance level desired for simulation


MSE <- function(recon, target, calibInt, MC=NULL, alpha=.01){
  
  #check that series are of equal length
  reconLength <- length(recon)
  if (sum(!is.na(recon)) != sum(!is.na(target))){
    stop("The two series must be of the same length!\n")
  }
  
  verifInt <- (1:length(target))[((1:length(target)) %in% calibInt[1]:calibInt[2]) == FALSE]
  
  #Show user the inputs they provided:
  cat("Length of common interval: ", sum(!is.na(recon)), "\n")
  cat("Calibration interval: ", calibInt[1], " - ", calibInt[2], "\n")
  cat("Verification interval: ", verifInt[1], " - ", verifInt[length(verifInt)], "\n\n")
  
  
  #Initiate variables and run RE/CE
  sumErr <- 0
  sumErrAvgV <- 0
  sumErrAvgC <- 0
  targetAvgC <- mean(target[calibInt[1]:calibInt[2]])
  targetAvgI <- mean(target)
  
  for (i in 1:length(recon)){
    #MSE calcs
    err <- (recon[i] - target[i])^2
    sumErr <- sumErr + err
    
    #RE calcs
    errAvgC <- (targetAvgC - target[i])^2
    sumErrAvgC <- sumErrAvgC + errAvgC
    
    #CE calcs
    errAvgV <- (targetAvgI - target[i])^2
    sumErrAvgV <- sumErrAvgV + errAvgV
  }
  
  mse <- sumErr/length(recon)
  mseC <- sumErrAvgC/length(recon)
  RE <- 1-(mse/mseC)
  mseV <- sumErrAvgV/length(recon)
  CE <- 1 - (mse/mseV)
  
  Rsquared <- cor(recon, target)^2
  
  
  if (!is.null(MC)){
    
    cat("Monte Carlo iterations for simulated significance level: ", MC, "\n")
    cat("Significance Level: ", alpha, "\n\n")
    
    #Repeat run of RE/CE to attain a significance Threshold
    if (require(astrochron) == F){
      install.packages("astrochron")
      library(astrochron)
    }
    
    RE2 <- rep(NA, MC)
    CE2 <- rep(NA, MC)
    RsquaredSim <- rep(NA, MC)
    
    #Generate surrogate time series of recon
    surrogate <- surrogates(target, nsim = MC, verbose = F, genplot = F)
    for (j in 1:MC){
      surrogate1 <- surrogate[,j]
      
      sumErr2 <- 0
      sumErrAvg2 <- 0
      sumErrAvgC2 <- 0
      sumErrAvgV2 <- 0
      
      for (i in 1:length(target)){
        err2 <- (surrogate1[i] - target[i])^2
        sumErr2 <- sumErr2 + err2
        
        errAvgC2 <- (targetAvgC - target[i])^2
        sumErrAvgC2 <- sumErrAvgC + errAvgC2
        
        errAvgV2 <- (targetAvgI - target[i])^2
        sumErrAvgV2 <- sumErrAvgV2 + errAvgV2
      }
      mse2 <- sumErr2/length(surrogate1)
      mseC2 <- sumErrAvgC2/length(surrogate1)
      RE2[j] <- 1-(mse2/mseC2)
      mseV2 <- sumErrAvgV2/length(surrogate1)
      CE2[j] <- 1-(mse2/mseV2)
      
      RsquaredSim[j] <- cor(surrogate1, target)^2
    }
    place <- MC - alpha*MC
    SLev1 <- sort(RE2)[place]
    SLev2 <- sort(CE2)[place]
    RsquaredSLev <- sort(RsquaredSim)[place]
    
    returns <- list("MSE" = mse, "RE" = RE, "SigLevRE" = SLev1, "CE" = CE, "SigLevCE" = SLev2,
                    "R2" = Rsquared, "SigLevR2" = RsquaredSLev, "MSEv" = mseV)
  }else{
    returns <- list("MSE" = mse, "RE" = RE, "CE" = CE, "R2" = Rsquared, "MSEv" = mseV)
  }
  
  return(returns)
}
