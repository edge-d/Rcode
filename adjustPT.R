#Take the PT rwl and original rwl files to adjust the relative positions of the PT series relative to the other series
#Find the relative distance of each measurement to the respective age-aligned average
##For each series, calculate the relative distance from the mean by averaging each individual offset
#repeat the previous two steps for the PT series
#Correct the placement of the PT series using the ratio
#Continue checking and rerunning the correction until the maximum error of an individual series is less than 10%
#Return the corrected PT rwl

library(rowr)

adjustPT <- function(rwlPT, rwlOrig, poFile){
  if(dim(rwlPT)[1] != dim(rwlOrig)[1]){
    stop("The two rwl files must be of the same dimensions")
  }
  if(dim(rwlPT)[2] != dim(rwlOrig)[2]){
    stop("The two rwl files must be of the same dimensions")
  }
  
  deltaPos <- 1
  
  while(max(abs(deltaPos)) > .1){
    #Age Align the original and PT data
    ageAligned <- NULL
    for (i in 1:length(colnames(rwlOrig))){
      naRemoved <- array(na.omit(rwlOrig[,i]))
      naFill <- rep(NA, poFile[i,2])
      newSeries <- c(naFill, naRemoved)
      ageAligned <- cbind.fill(ageAligned, newSeries, fill = NA)
    }
    ageAligned[,1] <- NULL
    colnames(ageAligned) <- colnames(rwlOrig)
    ###
    ageAlignedPT <- NULL
    for (i in 1:length(colnames(rwlPT))){
      naRemoved <- array(na.omit(rwlPT[,i]))
      naFill <- rep(NA, poFile[i,2])
      newSeries <- c(naFill, naRemoved)
      ageAlignedPT <- cbind.fill(ageAlignedPT, newSeries, fill = NA)
    }
    ageAlignedPT[,1] <- NULL
    colnames(ageAlignedPT) <- colnames(rwlPT)
    
    aaMean <- rowMeans(ageAligned, na.rm = TRUE)
    PTaaMean <- rowMeans(ageAlignedPT, na.rm = TRUE)
    ageAlignedRatio <- ageAligned
    PTageAlignedRatio <- ageAlignedPT
    for (i in 1:dim(ageAligned)[1]){
      #for i each row
      for (j in 1:dim(ageAligned)[2]){
        #for j each column
        if(!is.null(ageAligned[i,j])){
          #if there's a value at this point - row i, column j
          #find the ratio by dividing the age-aligned ring-width value by the row mean
          ageAlignedRatio[i,j] <- ageAligned[i,j]/aaMean[i]
          PTageAlignedRatio[i,j] <- ageAlignedPT[i,j]/PTaaMean[i]
        }
      }
      #relativePos[i] <- mean(ageAlignedRatio[i,], na.rm = TRUE)
    }
    relativePos <- colMeans(ageAlignedRatio, na.rm = TRUE)
    PTPos <- colMeans(PTageAlignedRatio, na.rm = TRUE)
    deltaPos <- relativePos - PTPos
    deltaPos <- deltaPos * mean(colMeans(rwlPT, na.rm = TRUE))
    
    for (i in 1:dim(rwlPT)[2]){
      rwlPT[,i] <- rwlPT[,i] + deltaPos[i]
    }
  }
  return(rwlPT)
}


