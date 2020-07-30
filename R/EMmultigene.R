## EMmultigene.R
## Functions for running expectation maximisation on data from multiple genes.
## Must use windows of particular size and calculate % methylation of
## methylatable sites in those windows.



#' Select reads from multiple genes
#' @param dataMatrix  A matrix containing the sample. dataMatrix [i,j] is
#' the methylation value or bincount value of sample i at position j.
#' @param minReads The minimum number of reads a matrix needs to have to be
#' sampled from.
#' @param addToReadName A string that will be added to the read names (row names)
#' of the matrix
#' @param preferBest If TRUE will sample reads in inverse proportion to the
#' fraction of NAs in the read (takes more informative reads preferentially)
selectReadsFromMatrix<-function(dataMatrix,minReads=50,addToReadName="",
                             preferBest=T){
  if(nrow(dataMatrix)>minReads){
    # append addtoReadName to row names:
    row.names(dataMatrix)<-paste(addToReadName,row.names(dataMatrix),sep="__")
    # sample in inverse proportion to number of NAs (preferably get the good
    # reads)
    if(preferBest) {
      fractionNA<-rowSums(is.na(dataMatrix))/ncol(dataMatrix)
      samplingWeight<-(1-fractionNA)/sum(1-fractionNA)
    } else {
      samplingWeight<-1/nrow(dataMatrix)
    }
    idx<-sample(1:nrow(dataMatrix),size=minReads,replace=F,prob=samplingWeight)
    return(dataMatrix[idx,])
  } else {
    return(NULL)
  }
}

