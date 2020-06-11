makeIDs <- function(rwi){
  source("E:/Lab Backup/R/is_odd.R")
  newIDS <- data.frame(matrix(ncol = 2, nrow = dim(rwi)[2]))
  colnames(newIDS) <- c("tree", "core")
  #Get rwi stats
  for (i in 1:dim(rwi)[2]){
    #step through i=1 and i=2 setting the "tree" ID to 001 and the "core" ID to 1 and 2
    newIDS[i,"tree"] <- ceiling(i/2)
    if(is.odd(i)){
      newIDS[i,"core"] <- 1
    }else(
      newIDS[i,"core"] <- 2
    )
  }
  return(newIDS)
}
  