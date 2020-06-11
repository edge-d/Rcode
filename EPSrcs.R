#Function returns running rwi statistics from the dplR function "rwi.stats.running"
#EPS statistic is adjusted for increased variance of RCS based on methods of Melvin and Briffa, 2014
#requires package dplR and function "makeIDs"
#Warning, makeIDs assumes 2 cores per tree
#inputs: rwl = an rwl file as produced by dplR's read.rwl()
# rwi = an rwi file as produced by dplR's detrend() and corresponding to the rwl file
# ... = additional inputs to dplR's rwi.stats.running

EPSrcs <- function(rwl, rwi,  ...){
  library(dplR)
  source("E:/Lab Backup/R/makeIDs.R")
  
  #Get original rwi stats
  RCSrwiStats <- rwi.stats.running(rwi, ids = makeIDs(rwi), window.length = ...)
  
  repeats1 <- length(RCSrwiStats$start.year)
  EPSadj <- rep(NA, repeats1)
  for (i in 1:repeats1){
    #calculate Nadj for RCS
    #Nadj = n x (50-year spline variance / RCS variance)
    splineRWI <- detrend(rwl, method = "Spline", nyrs = 50)
    splineRWIsub <- splineRWI[as.numeric(rownames(splineRWI)) %in%  RCSrwiStats$start.year[i]:RCSrwiStats$end.year[i],]
    splineSD <- mean(apply(splineRWIsub, 1, function(x) sd(x, na.rm = TRUE)))
    RCSrwisub <- rwi[as.numeric(rownames(rwi)) %in% RCSrwiStats$start.year[i]:RCSrwiStats$end.year[i],]
    RCSsd <- mean(apply(RCSrwisub, 1, function(x) sd(x, na.rm = TRUE)))
    Nadj <- RCSrwiStats$n[i] * splineSD/RCSsd
    EPSadj[i] <- round((Nadj * RCSrwiStats$rbar.eff[i])/((Nadj * RCSrwiStats$rbar.eff[i]) + (1-RCSrwiStats$rbar.eff[i])), 3)
  }
  RCSrwiStats <- cbind(RCSrwiStats, EPSadj)
  return(RCSrwiStats)
}