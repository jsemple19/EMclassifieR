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
#' fraction of NAs in the read (takes more informative reads preferentially).
#' This is calculated using the fourth power of the fraction of NAs in order to
#' bias more strongly for better matrices.
#' @return Returns a dataMatrix with minReads number of rows, subsampled from
#' the original dataMatrix.
#' @export
selectReadsFromMatrix<-function(dataMatrix, minReads=50, addToReadName="",
                                preferBest=T){
  if(nrow(dataMatrix)>minReads){
    # append addtoReadName to row names:
    row.names(dataMatrix)<-paste(addToReadName,row.names(dataMatrix),sep="__")
    # sample in inverse proportion to number of NAs (preferably get the good
    # reads)
    if(preferBest) {
      fractionNA<-rowSums(is.na(dataMatrix))/ncol(dataMatrix)
      samplingWeight<-(1-fractionNA)^4/sum((1-fractionNA)^4)
    } else {
      samplingWeight<-1/nrow(dataMatrix)
    }
    idx<-sample(1:nrow(dataMatrix),size=minReads,replace=F,prob=samplingWeight)
    return(dataMatrix[idx,])
  } else {
    return(NULL)
  }
}


#' Get full matrix of all positions in a window even if no C
#'
#' In order to compare different promoters, we need to create a padded
#' matrix with NAs in positions in between Cs.
#'
#' @param dataMatrix A matrix of methylation values (reads x position) with
#' positions being relative coordinates to an anchor point
#' @param winSize The size (in bp) of the window containing the matrices
#' @param anchorPoint One of "middle" or "start": the position from which numbering starts for the positions in the dataMatrix
#' @return A padded methylation matrix
#' @export
getFullMatrix<-function(dataMatrix,winSize=500,anchorPoint="middle") {
  fullMat<-matrix(data=NaN,nrow=dim(dataMatrix)[1],ncol=winSize)
  Cpos<-colnames(dataMatrix)
  if(anchorPoint=="middle") {
    withinRange<- -winSize/2<=as.numeric(Cpos) & winSize/2>=as.numeric(Cpos)
    colnames(fullMat)<-c(as.character(seq(-winSize/2,-1)),
                         as.character(seq(1,winSize/2)))
  } else if(anchorPoint=="start") {
    withinRange<- 0<=as.numeric(Cpos) & winSize>=as.numeric(Cpos)
    colnames(fullMat)<-as.character(seq(1,winSize))
  } else {
    print("anchor point should be one of 'start' or 'middle'")
  }
  fullMat[,Cpos[withinRange]]<-dataMatrix[,withinRange]
  rownames(fullMat)<-rownames(dataMatrix)
  return(fullMat)
}



#' Make window matrix
#'
#' @param fullMatrix  A matrix containing methylation values for a full
#' region of the genome using coordinates relative to an anchor point, and
#' includes non-Cs that are marked by NaN. dataMatrix [i,j] is the methylation
#' value of sample i at position j.
#' @param binSize Size of window (in bp) in which to average results.
#' @param stepSize size of step by which to slide the window forward.
#' @return Matrix with average fraction methylation within windows. Normalised
#' by the number of informative postions. NaN indicates a window with no
#' methylatable positions. NA indicates a window with no informative methylation
#' call.
#' @export
prepareWindows <- function(fullMatrix, binSize=20, stepSize=1) {
  print(paste("Preparing windows for window size", binSize))
  binStarts<-seq(1, ncol(fullMatrix)-binSize+1, by=stepSize)
  binEnds<-seq(binSize, ncol(fullMatrix), by=stepSize)
  winMatrix<-matrix(data=NaN, nrow=nrow(fullMatrix),
                    ncol=ncol(fullMatrix)-binSize+1)
  startVals<-as.numeric(colnames(fullMatrix)[binStarts])
  endVals<-as.numeric(colnames(fullMatrix)[binEnds])
  midPoints<-floor(rowMeans(cbind(startVals,endVals)))
  colnames(winMatrix)<-as.character(midPoints)
  rownames(winMatrix)<-rownames(fullMatrix)
  #Calculate values for each bin
  for(i in 1:length(binStarts)){
    noSites<-rowSums(is.nan(fullMatrix[,i:(i+binSize-1)]))==binSize
    informative<-rowSums(!is.na(fullMatrix[,i:(i+binSize-1)]))
    methylated<-rowSums(fullMatrix[,i:(i+binSize-1)],na.rm=T)
    winMatrix[,i]<-ifelse(noSites,NaN,
                          ifelse(!informative, NA, methylated / informative))
    #winMatrix[,i]<-ifelse(!informative, NA, round(methylated/informative,2))
  }
  return(winMatrix)
}



### to convert scale to -1 to 1, do 2x-1, then code NA and NaN as 0, or NaN as 0,
### and NA as randomly assigned.

#' Recode matrix
#'
#' Recode values in methylation matrix from a scale of 0 (no meth) to 1 (meth)
#' to a scale of -1 (no meth) to 1 (meth). NaN are set to 0 and NAs can
#' be either 0 or randomly assigned to -0.5 or 0.5
#' @param dataMatrix  A matrix containing methylation values for a
#' region of the genome using coordinates relative to an anchor point.
#' dataMatrix[i,j] is the methylation value of sample i at position j.
#' 0=not methylated, 1=methylated, NA=no methylation information at this
#' position, NaN=non-methylatable positon.
#' @param randomiseNAs if TRUE NAs randomly receive values between -0.5 and 0.5.
#' If FALSE, NAs are encoded as 0
#' @return Matrix with methylation values coded on a -1 to 1 scale
#' @export
rescale_minus1To1<-function(dataMatrix,randomiseNAs=F){
  dataMatrix[!is.na(dataMatrix)]<-dataMatrix[!is.na(dataMatrix)]*2-1
  dataMatrix[is.nan(dataMatrix)]<-0
  if(randomiseNAs) {
    dataMatrix[is.na(dataMatrix)]<-sample(seq(-0.5,0.5,by=0.01),
                                          sum(is.na(dataMatrix)),
                                          replace=T)
  } else {
    dataMatrix[is.na(dataMatrix)]<-0
  }
  return(dataMatrix)
}



#' #' Recode matrix
#' #'
#' #' Recode values in methylation matrix from a scale of 0 (no meth) to 1 (meth)
#' #' to a scale of -1 (no meth) to 1 (meth). NaN are set to 0 and NAs can
#' #' be either 0 or randomly assigned to -0.5 or 0.5
#' #' @param dataMatrix  A matrix containing methylation values for a
#' #' region of the genome using coordinates relative to an anchor point.
#' #' dataMatrix[i,j] is the methylation value of sample i at position j.
#' #' 0=not methylated, 1=methylated, NA=no methylation information at this
#' #' position.
#' #' @return Matrix with methylation values coded on a -1 to 1 scale
#' #' @export
#' recodeMatrixAsNumeric<-function(dataMatrix){
#'   dataMatrix[is.na(dataMatrix)]<-0.5
#'   return(dataMatrix)
#' }

### to convert scale to -1 to 1, do 2x-1, then code NA and NaN as 0, or NaN as 0,
### and NA as randomly assigned.

#' Recode matrix
#'
#' Recode values in methylation matrix from a scale of 0 (no meth) to 1 (meth)
#' to a scale of -1 (no meth) to 1 (meth). NaN are set to 0 and NAs can
#' be either 0 or randomly assigned to -0.5 or 0.5
#' @param dataMatrix  A matrix containing methylation values for a
#' region of the genome using coordinates relative to an anchor point.
#' dataMatrix[i,j] is the methylation value of sample i at position j.
#' -1=not methylated, 1=methylated. NAs and NaNs have been assigned the value 0
#' @return Matrix with methylation values coded on a -1 to 1 scale
#' @export
rescale_0To1<-function(dataMatrix){
  dataMatrix<-(dataMatrix+1)/2
  return(dataMatrix)
}