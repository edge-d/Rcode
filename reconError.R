#Bootstrapping simulation to estimate prediction intervals of the reconstructed SST by: 
#1) we used the periods 1941:1971 and 1972:2002 as independent calibration and verification intervals, respectively, 
#2) an SST reconstruction model was developed by linear regression in the calibration interval, 
#3) at each possible sample depth 1:n, where n is the sample depth of the year in the verification period with the fewest replicates, 
##a) for each iteration i in 1:1000 the detrended indices in the verification period were resampled with replacement (the resampling 
##disregards the continuity of indices such that year 1 and year 2 of a single resampling may be represented by different individuals),
###i) resampled indices were combined by robust bi-weight mean, 
###ii) an SST reconstruction was developed for the verification interval based on the regression coefficients from step 2, 
###iii) each reconstruction was aligned with the instrumental SST and the difference (absolute value) was calculated between the 
##instrumental SST and reconstructed value at each year,  
###iv) the ascending-order 90th percentile value was stored, 
##b) the median of all 90th percentile error values at a given sample depth is stored (error values ± reconstructed values represent 
##the 90% confidence envelope), 
#4) steps 2 and 3 were repeated with the calibration and verification intervals reversed, 
#5) the mean of the 90th percentile values at each sample depth from the two calibration-verification iterations was calculated
#6) the bootstrapped 90th percentile error was applied to the SST reconstruction as a prediction interval by direct addition to 
#(and subtraction from) the reconstructed SSTs prior to the instrumental period ().


################################################################################################

#Inputs Required:

  #Period of overlap
    #start and end years of chronology-target overlap interval as in c(1941,2002)
  #chronology rwi
    #rwi as built by dplR
  #reconstruction target
    #two-column data frame with "Year" as first column name

#Optional inputs:
  #numIts
    #Number of iterations to resample at a given sample depth, defaults to 1000
  #LogTrans
    #Should the chronology be log transformed before being compared with the target (T/F)
  #verifInt
      #Length of the verification interval in years, defaults to 10
  #CI
    #The desired prediciction Interval as decimal (ex. 95%, input as 0.95), defaults to 95%


################################################################################################

#Reguires packages: dplR
#requires functions: qErr(), cbind.NA()

reconError <- function(rwi, targetDat, pOverlap, numIts = 1000, CI = .95, logTrans = FALSE, verifInt = 10){
  #Get packages
  library(dplR)
  #Import functions
  source("E:/Lab Backup/R/qErr.R")
  source("E:/Lab Backup/R/cbind_NA.R")

  #Truncate rwl to overlap interval
  rwi <- rwi[as.numeric(rownames(rwi)) %in% pOverlap[1]:pOverlap[2],]
  ################################################################################################
  #Split period into calib and verif
  endP1 <- pOverlap[2] - verifInt
  startP2 <- endP1 + 1
  P1 <- c(pOverlap[1]:endP1)
  P2 <- c(startP2:pOverlap[2])
  
  for (Period in 1:2){
      
    ################################################################################################
    #Get chronology for calib interval
    rwi1 <- rwi[as.numeric(rownames(rwi)) %in% P1,]
    P1chron <- chron(rwi1)
    #Run regression for SST model
    df2 <- cbind.data.frame(P1chron$xxxstd, targetDat$SST[targetDat$Year %in% P1])
    colnames(df2) <- c("chron", "target")
    if(logTrans == TRUE){
      df2$chron <- log(df2$chron+10)
    }
    #Get recon ceofs
    coefs <- summary(lm(target~chron, data = df2))
    coefs <- coefs$coefficients
    
    ################################################################################################
    #Truncate rwi to series with coverage for full verif interval
    rwi2 <- rwi[as.numeric(rownames(rwi)) %in% P2,]

    #Get target for verif period
    targetVerif <- targetDat[targetDat$Year %in% P2,2]
    
    #Total number of samples that cover the full verification interval
    sampleMax <- min(apply(rwi2, 1, function(x) sum(!is.na(x))))
    #Initialize data frame to collect RMSE for each sample depth for all iterations
    RMSEfull <- data.frame(matrix(nrow = sampleMax, ncol = numIts, data = NA))
    rownames(RMSEfull) <- as.character(1:sampleMax)
    #Initialize sampling vector
    newSamples <- data.frame()
    
    for (i in 1:sampleMax){
      for (j in 1:numIts){
        #sample from verif interval based on sample depth
        
        for (k in 1:dim(rwi2)[1]){
          valuesAtYeark <- rwi2[k,!is.na(rwi2[k,])]
          newSamples[k,1:i] <- sample(valuesAtYeark, i, replace = TRUE)
        }
        #Build chrono
        if(i != 1){
          sampleChron <- chron(newSamples)
          if(logTrans == TRUE){
            sampleChron$xxxstd <- log(sampleChron$xxxstd+10)
          }
          #Build verif period recon
          sampleRecon <- (sampleChron$xxxstd * coefs[2,1]) + coefs[1,1]
        }else{#special case for sample depth = 1
          if(logTrans == TRUE){
            newSamples <- log(newSamples+10)
          }
          sampleRecon <- (newSamples * coefs[2,1]) + coefs[1,1]
        }
        
        #Calculate RMSE in verif interval
        RMSEfull[i,j] <- qErr(sampleRecon, targetVerif, CI)
      }
    }
    
    if(Period ==2){
      newRMSEdat <- data.frame(apply(RMSEfull, 1, function(x) quantile(x, c(.5))))
      RMSEdat <- cbind.NA(RMSEdat, newRMSEdat)
      colnames(RMSEdat) <- c("LateVerif", "EarlyVerif")
      rownames(RMSEdat) <- as.character(1:dim(RMSEdat)[1])
    }else if (Period ==1){
      RMSEdat <- data.frame(apply(RMSEfull, 1, function(x) quantile(x, c(.5))))
      rownames(RMSEdat) <- as.character(1:sampleMax)
      colnames(RMSEdat) <- "LateVerif"
      P3 <- P1
      P1 <- P2
      P2 <- P3
    }
    
  }

  
  return(RMSEdat)

}
