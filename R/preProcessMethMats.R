## functions for preprocessing methylation matrices


#' Remove rows containing NAs
#'
#' @param mat A matrix of numbers and NAs.
#' @return A matrix without the rows containing NAs
#' @examples
#' removeAllNArows(matrix(c(1,2,3,NA,4,5),nrow=3,byrow=T))
#' @export
removeAllNArows<-function(mat) {
  NAidx<-rowSums(is.na(mat))>0
  newMat<-mat[!NAidx,]
  return(newMat)
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
#' padMissingColWithNAs(matrix(c(1:4),nrow=2,byrow=T,dimnames=list(c("a","b"),c(-1,2))))
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
