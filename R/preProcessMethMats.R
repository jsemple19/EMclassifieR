## functions for preprocessing methylation matrices


#' Remove rows containing NAs
#'
#' Function to remove all rows from a matrix where more than a certain
#' fraction of the positions in that row are NAs. To remove rows with any
#' NAs at all, set maxNAfraction=0
#' @param dataMatrix A matrix of numbers and NAs.
#' @param maxNAfraction Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)
#' @param removeAll0 Remove reads that only have 0 or NA values, i.e. were not methylated at all (default=F)
#' @return A matrix without the rows where the fraction of NA positions is above
#' the threshold
#' @examples
#' removeNArows(matrix(c(1,2,3,NA,4,5),nrow=3,byrow=TRUE))
#' @export
removeNArows<-function(dataMatrix, maxNAfraction=0.2, removeAll0=F) {
  NAidx<-rowSums(is.na(dataMatrix))/dim(dataMatrix)[2] > maxNAfraction
  newMat<-dataMatrix[!NAidx,]
  if(removeAll0){
    idx0orNA<-rowSums(dataMatrix==0 | is.na(dataMatrix))==dim(dataMatrix)[2]
    newMat<-dataMatrix[!idx0orNA,]
  }
  return(newMat)
}


#' Recode matrix
#'
#' Recode values in methylation matrix from a scale of 0 (no meth) to 1 (meth)
#' to a scale of -1 (no meth) to 1 (meth). NaN are set to 0 and NAs can
#' be either 0 or randomly assigned to -0.5 or 0.5
#' @param dataMatrix  A matrix containing methylation values for a
#' region of the genome using coordinates relative to an anchor point.
#' dataMatrix[i,j] is the methylation value of sample i at position j.
#' 0=not methylated, 1=methylated, NA=no methylation information at this
#' position.
#' @return Matrix with methylation values coded on a -1 to 1 scale
#' @export
recodeMatrixAsNumeric<-function(dataMatrix){
  dataMatrix[is.na(dataMatrix)]<-0.5
  return(dataMatrix)
}




#' Add padding between informative positions.
#'
#' In order to combine matrices from different genes, the simple
#' read x position matrix must be converted to a full matrix with
#' a column for every position within the region, even if no CpG or GpC
#' is found there in a particular gene. Any non-informative position gets
#' the value "NA".
#'
#' @param mat A matrix of numbers and NAs.
#' @param colRange a vector of two numbers for start and end positions used in the matrix (default=c(-250,250))
#' @return A matrix with missing columns padded with NAs
#' @examples
#' m<-matrix(c(1:4),nrow=2,byrow=TRUE,dimnames=list(c("a","b"),c(-1,2)))
#' padMissingColWithNAs(m,colRange=c(-2,2))
#' @export
padMissingColWithNAs<-function(mat,colRange=c(-250,250)) {
  notNA<-colnames(mat)
  newColNumber=colRange[2]-colRange[1]
  fullMat<-matrix(rep(NA,nrow(mat)*newColNumber),ncol=newColNumber)
  colnames(fullMat)<-remove0(colRange[1]:colRange[2])
  fullMat[,colnames(mat)]<-mat
  rownames(fullMat)<-rownames(mat)
  return(fullMat)
}



#' Remove the number 0 from a vector
#'
#' @param vec A vector of numbers
#' @return The same vector without the number 0
#' @examples
#' remove0(c(-2:2))
#' @export
remove0<-function(vec) {
  i<-match(0,vec)
  if (!is.na(i)) {
    vec<-vec[-i]
  }
  return(vec)
}
