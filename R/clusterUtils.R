# functions to test clustering quality

#' Calculate the sum of square distances between rows of a matrix
#'
#' @param dataMatrix A matrix for which to calcluate the sum of the sum of square distances between its rows.
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return Numeric. The sum of the square distances between the rows of the matrix.
#' @export
SumSqMatrixRows<-function(dataMatrix,distMetric=list(name="euclidean")) {
  #sum(stats::dist(dataMatrix)^2)
  sum(getDistMatrix(dataMatrix,distMetric)^2)
}


#' Calculate within-cluster sum of squares
#'
#' @param dataMatrix A matrix for which you wish to calculate the within cluster sum of squares distance between rows
#' @param classes A vector of classes designating which rows belong to which class
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return A number indicating the sum of all within-cluster variation
#' @export
withinClusterSS<-function(dataMatrix,classes,distMetric=list(name="euclidean")){
  SS<-0
  for (class in unique(classes)) {
    tmpMat<-dataMatrix[classes==class,]
    if (!is.null(dim(tmpMat))) {
      SS<-SS+SumSqMatrixRows(tmpMat,distMetric)
    }
  }
  return(SS)
}



#' Randomise the values in each row of a matrix
#'
#' @param dataMatrix A matrix to be randomised
#' @return A matrix where each row has been randomised
#' @export
randomiseMatrixRows<-function(dataMatrix){
  t(apply(dataMatrix,1,function(r){
    sample(r,replace=F)
  }))
}


#' Estimate within cluster sum of squares
#'
#' Get estimated total within cluster sum of squares by clustering random matrices
#' @param dataMatrix Matrix to be randomised and clustered
#' @param k_range A vector indicating different numbers of classes to learn
#' @param maxB  The maximum number of randomisations to perform
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param nThreads Number of threads to use for generating background distribution (default is 1)
#' @param setSeed Logical value to determine if seed should be set for randomisation (default is FALSE)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return A data frame with the average of the total within class sum of squares for multiple randomised matrices and different numbers of classes
#' @export
clusterRandomMatrices<-function(dataMatrix, k_range=2:8, maxB=100,
                                convergenceError=1e-6, maxIterations=100,
                                nThreads=1, setSeed=F,
                                distMetric=list(name="euclidean")){
  totalWSS<-data.frame(numClasses=k_range, meanWSS=0, sumSq=0, sdWSS=NA)
  #clst<-parallel::makeCluster(nThreads)
  #doParallel::registerDoParallel(clst)
  #if(setSeed) { doRNG::registerDoRNG(123) }
  if(isTRUE(setSeed)){
    set.seed(200413)
  }
  for (numClasses in k_range){
    print(paste0("numClasses: ",numClasses))
    nc<-which(totalWSS$numClasses==numClasses)
   # allwss<-foreach::foreach(1:maxB, .combine=c,.export=c("randomiseMatrixRows",
   #                                    "runEMrepeats_withinSS")) %dopar%{
      for (b in 1:maxB){

      randMat<-randomiseMatrixRows(dataMatrix)
      wss<-runEMrepeats_withinSS(randMat, numClasses, convergenceError,
                                maxIterations, EMrepeats=1,distMetric)
      totalWSS[nc,"meanWSS"] <- totalWSS[nc,"meanWSS"]+wss/maxB
      totalWSS[nc,"sumSq"] <- totalWSS[nc,"sumSq"]+(wss^2)/maxB
    }
    #totalWSS[nc,"meanWSS"] <- mean(allwss)
    #totalWSS[nc,"sumSq"]<-sum(allwss^2/maxB)
    totalWSS[nc,"sdWSS"] <- sqrt(totalWSS[nc,"sumSq"]-totalWSS[nc,"meanWSS"]^2)
  }
  #parallel::stopCluster(clst)
  return(totalWSS)
}



#' Test if matrix is valid
#'
#' Test if matrix contains values that are only within the range specified by valueRange
#' @param dataMatrix Matrix the values of which will be tested
#' @param valueRange Vector of min and max permissible values in the matrix (default=c(0,1))
#' @param NAsValid Boolean TRUE/FALSE indicating if NAs are considered valid (default=FALSE)
#' @return Boolean TRUE/FALSE indicating if matrix contains only valid values
#' @export
isMatrixValid<-function(dataMatrix,valueRange=c(0,1),NAsValid=FALSE){
  isValid=TRUE
  if(!is.matrix(dataMatrix)) {
    print("Not a proper matrix")
    isValid=FALSE
  }
  if(is.matrix(dataMatrix) & any(dim(dataMatrix)==0)) {
    print("Empty dimension")
    isValid=FALSE
  }
  if(!NAsValid & sum(is.na(dataMatrix))>0) {
    print("Matrix contains NA values")
    isValid=FALSE
  }
  NAidx<-is.na(dataMatrix)
  if(sum(dataMatrix[!NAidx]<valueRange[1])>0) {
    print(paste("Matrix contains values less than",valueRange[1]))
    isValid=FALSE
  }
  if(sum(dataMatrix[!NAidx]>valueRange[2])>0) {
    print(paste("Matrix contains values greater than",valueRange[2]))
    isValid=FALSE
  }
  return(isValid)
}


#' Calculate Euclidean distance between two vectors
#'
#' @param r1 A numeric vector
#' @param r2 Second numeric vector of same length as r1
#' @return Euclidean distance between two vectors
#' @export
euclidean<-function(r1,r2){
  stopifnot(length(r1)==length(r2))
  sqrt(sum((r1-r2)^2))
}


#' Calculate Euclidean distance between two vectors with sliding windows
#'
#' To increase dynamic range and spatial information in comparin two binary
#' vectors, a sliding window is used. The euclidean distance is computed within
#' each window, and then moved by a step of 1.
#' @param r1 A numeric vector
#' @param r2 Second numeric vector of same length as r1
#' @param winSize Sliding window size (number of values to combine)
#' @return Euclidean distance between two vectors
#' @export
euclidWin<-function(r1,r2,winSize=3){
  stopifnot(length(r1)==length(r2), length(r1)>=winSize)
  distSum<-0
  for(i in 1:(length(r1)-winSize+1)){
    distSum<-distSum+euclidean(r1[i:(i+winSize-1)],r2[i:(i+winSize-1)])
    print(distSum)
  }
  return(distSum)
}




#' Calculate Euclidean distance between all rows of a matrix with sliding window
#'
#' To increase dynamic range and spatial information in comparing a matrix of
#' vectors, a sliding window is used. The euclidean distance is computed within
#' each window, and then moved by a step of 1.
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' @param winSize Sliding window size (number of values to combine)
#' @param nThreads Number of threads to use for generating background distribution (default is 1)
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
euclidWinDist<-function(binMat,winSize=3,nThreads=1){
  i <- j <- NULL
  numRows<-nrow(binMat)
  distMat<-matrix(rep(NA,numRows^2),nrow=numRows,ncol=numRows)
  clst<-parallel::makeCluster(nThreads)
  doParallel::registerDoParallel(clst)
  distList<-foreach::foreach(i=2:numRows,.combine=c,.export=c("euclidWin",
                                                              "euclidean")) %:%
    foreach::foreach(j=1:(i-1),.combine=c) %dopar%{
  #for(i in 1:numRows){
    #for(j in 1:i){
      #distMat[i,j]<-euclidWin(binMat[i,],binMat[j,])
      euclidWin(binMat[i,],binMat[j,],winSize=winSize)

    #}
    }
  parallel::stopCluster(clst)
  distMat[upper.tri(distMat)]<-t(distList)
  distMat<-t(distMat)
  rownames(distMat)<-rownames(binMat)
  colnames(distMat)<-rownames(binMat)
  return(stats::as.dist(distMat))
}




#' A generic function for implementing various distance metrics
#'
#' To avoid writing too much code when i want to try a different distance
#' function i want a generic function that will be set up appropriately
#' with the right parameters (not specified before)
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' @param distMetric List continaing the name of the metric and the values of
#' any parameters the metric nromally takes.
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
getDistMatrix<-function(binMat,distMetric=list(name="euclidean")){
  switch(distMetric$name,
         euclidean=stats::dist(binMat),
         euclidWinDist=euclidWinDist(binMat,winSize=distMetric$winSize,
                                     nThreads=distMetric$nThreads))
}
