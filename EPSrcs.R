#Function returns running rwi statistics from the dplR function "rwi.stats.running"
#EPS statistic is adjusted for increased variance of RCS based on methods of Melvin and Briffa, 2014
#https://doi.org/10.1016/j.dendro.2014.07.008
#requires package dplR and function "makeIDs"
#Warning, makeIDs assumes 2 cores per tree
#inputs: rwl = an rwl file as produced by dplR's read.rwl()
# rwi = an rwi file as produced by dplR's detrend() and corresponding to the rwl file
# winLen and WinOv correspond to dplR's rwi.stats.running inputs window.length and window.overlap
#spLen = length of spline for detrending to compare variance (original publcation 
#recommends 50-year spline, but the default is set to 2/3 series length) 

EPSrcs <- function(rwl, rwi, winLen = 30, winOv = 15, spLen = NULL, cores=1){
  library(dplR)
  source("D:/Lab Backup/R/makeIDs.R")
  
  #Get original rwi stats
  if(cores==1){
    RCSrwiStats <- rwi.stats.running(rwi, window.length = winLen, window.overlap = winOv)
  }else if (cores==2){
    RCSrwiStats <- rwi.stats.running(rwi, ids = makeIDs(rwi), window.length = winLen, window.overlap = winOv)
  }else{
    stop("Cores per tree must be either 1 or 2!")
  }
  
  
  repeats1 <- length(RCSrwiStats$start.year)
  EPSadj <- rep(NA, repeats1)
  Nadj <- rep(NA, repeats1)
  splineRWI <- dplR::detrend(rwl, method = "Spline", nyrs = spLen)
  for (i in 1:repeats1){
    #calculate Nadj for RCS
    #Nadj = n x (50-year spline variance / RCS variance)
    splineRWIsub <- splineRWI[as.numeric(rownames(splineRWI)) %in%  RCSrwiStats$start.year[i]:RCSrwiStats$end.year[i],]
    splineSD <- mean(apply(splineRWIsub, 1, function(x) sd(x, na.rm = TRUE)))
    RCSrwisub <- rwi[as.numeric(rownames(rwi)) %in% RCSrwiStats$start.year[i]:RCSrwiStats$end.year[i],]
    RCSsd <- mean(apply(RCSrwisub, 1, function(x) sd(x, na.rm = TRUE)))
    Nadj[i] <- RCSrwiStats$n[i] * splineSD/RCSsd
    EPSadj[i] <- round((Nadj[i] * RCSrwiStats$rbar.eff[i])/((Nadj[i] * RCSrwiStats$rbar.eff[i]) + (1-RCSrwiStats$rbar.eff[i])), 3)
  }
  RCSrwiStats <- cbind(RCSrwiStats, EPSadj, Nadj)
  return(RCSrwiStats)
}
