 # functions to test clustering quality

#' Calculate the sum of square distances between rows of a matrix
#'
#' @param dataMatrix A matrix for which to calcluate the sum of the sum of square distances between its rows.
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return Numeric. The sum of the square distances between the rows of the matrix.
#' @export
SumSqMatrixRows<-function(dataMatrix,distMetric=list(name="euclidean", rescale=F)) {
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
withinClusterSS<-function(dataMatrix,classes,distMetric=list(name="euclidean", rescale=F)){
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
                                distMetric=list(name="euclidean", rescale=F)){
  totalWSS<-data.frame(numClasses=k_range, meanWSS=0, sumSq=0, sdWSS=NA)
  clst<-parallel::makeCluster(nThreads)
  doParallel::registerDoParallel(clst)
  if(setSeed) { doRNG::registerDoRNG(123) }
  #if(isTRUE(setSeed)){
  #  set.seed(200413)
  #}
  for (numClasses in k_range){
    print(paste0("numClasses: ",numClasses))
    nc<-which(totalWSS$numClasses==numClasses)
    allwss<-foreach::foreach(1:maxB, .combine=c,.export=c("randomiseMatrixRows",
                                       "runEMrepeats_withinSS")) %dopar%{
    #  for (b in 1:maxB){
      browser()
      randMat<-randomiseMatrixRows(dataMatrix)
      wss<-runEMrepeats_withinSS(randMat, numClasses, convergenceError,
                                maxIterations, EMrepeats=1, distMetric)
    #  totalWSS[nc,"meanWSS"] <- totalWSS[nc,"meanWSS"]+wss/maxB
    #  totalWSS[nc,"sumSq"] <- totalWSS[nc,"sumSq"]+(wss^2)/maxB
    }
    totalWSS[nc,"meanWSS"] <- mean(allwss)
    totalWSS[nc,"sumSq"]<-sum(allwss^2/maxB)
    totalWSS[nc,"sdWSS"] <- sqrt(totalWSS[nc,"sumSq"]-totalWSS[nc,"meanWSS"]^2)
  }
  parallel::stopCluster(clst)
  return(totalWSS)
}



#' Test if matrix is valid
#'
#' Test if matrix contains values that are only within the range specified by valueRange
#' @param dataMatrix Matrix the values of which will be tested
#' @param valueRange Vector of min and max permissible values in the matrix
#' (default=c(0,1))
#' @param NAsValid Boolean TRUE/FALSE indicating if NAs are considered valid (default=FALSE)
#' @param checkColOrder Boolean TRUE/FALSE indicating if the order of column names should be checked. this assumes that these are numerical positions in ascending order (default=TRUE)
#' @return Boolean TRUE/FALSE indicating if matrix contains only valid values
#' @export
isMatrixValid<-function(dataMatrix, valueRange=c(0,1), NAsValid=FALSE,
                        checkColOrder=TRUE){
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
  if(checkColOrder){
    myColNames<-colnames(dataMatrix)
    if(any(myColNames!=as.character(sort(as.numeric(myColNames))))){
      print(paste("Matrix columns in wrong order"))
      isValid=FALSE
    }
  }
  return(isValid)
}



#' Calculate cosine distance between all rows of a matrix
#'
#'Using the cosine distance function form las package
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' @param valNA Value to give NAs in matrix
#' @param rescale Should matrix be rescaled from 0..1 range to -1..1 range?
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
cosineDist<-function(binMat, valNA=NULL, rescale=F){
  if(!is.null(valNA)){
    binMat[is.na(binMat)]<-valNA
  }
  if(rescale){
    binMat<-rescale_minus1To1(binMat)
  }
  print(paste0("Cosine clustering. NAs:", sum(is.na(binMat)),
               ". Matrix range:", min(binMat),"-", max(binMat)))
  cosDist<-stats::as.dist(1-lsa::cosine(t(binMat)))
  return(cosDist)
}



#' Calculate canberra distance between all rows of a matrix
#'
#'Using the canberra distance function from stats package
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' @param valNA Value to give NAs in matrix
#' @param rescale Should matrix be rescaled from 0..1 range to -1..1 range?
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
canberraDist<-function(binMat, valNA=NULL, rescale=F){
  if(!is.null(valNA)){
    binMat[is.na(binMat)]<-valNA
  }
  if(rescale){
    binMat<-rescale_minus1To1(binMat)
  }
  print(paste0("Canberra clustering. NAs:", sum(is.na(binMat)),
                ". Matrix range:", min(binMat),"-", max(binMat)))
  canDist<-stats::dist(binMat, method="canberra")
  return(canDist)
}


#' Calculate euclidean distance between all rows of a matrix
#'
#'Using the euclidean distance function from stats package
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' @param valNA Value to give NAs in matrix
#' @param rescale Should matrix be rescaled from 0..1 range to -1..1 range?
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
euclideanDist<-function(binMat, valNA=NULL, rescale=F){
  if(!is.null(valNA)){
    binMat[is.na(binMat)]<-valNA
  }
  if(rescale){
    binMat<-rescale_minus1To1(binMat)
  }
  print(paste0("Euclidean clustering. NAs:", sum(is.na(binMat)),
               ". Matrix range:", min(binMat),"-", max(binMat)))
  eucDist<-stats::dist(binMat, method="euclidean")
  return(eucDist)
}





#' A generic function for implementing various distance metrics
#'
#' To avoid writing too much code when i want to try a different distance
#' function i want a generic function that will be set up appropriately
#' with the right parameters (not specified before)
#' @param binMat A matrix of numbers for which you want to calculate the
#' distance between rows
#' #' @param valNA Value to give NAs in matrix
#' @param distMetric List continaing the name of the metric and the values of
#' any parameters the metric nromally takes.
#' @return A distance object (lower triangle) with the distances between all
#' rows of the input matrix
#' @export
getDistMatrix<-function(binMat, distMetric=list(name="euclidean",
                                                valNA=NULL, rescale=F)){
  tryCatch(
  {
    exists("distMetric$valNA", mode="numeric")
    valNA<-distMetric$valNA
  }, error=function(e){
    valNA<-NULL
    return(valNA)
  })

  switch(distMetric$name,
         euclidean=euclideanDist(binMat, valNA=valNA,
                                 rescale=distMetric$rescale),
         cosineDist=cosineDist(binMat, valNA=valNA,
                               rescale=distMetric$rescale),
         canberraDist=canberraDist(binMat, valNA=valNA,
                            rescale=distMetric$rescale))
}





#' Plot pca of data matrix
#'
#' @param dataMatrix Matrix for PCA with row names containing __classNumber at
#' the end
#' @param classes Vactor of classifications for all the rows
#' @param rescale Rescale matrix with 0-1 values to -1 to 1 values where -1 is
#' no methylation, NAs are 0 and 1 is methylation at that position.
#' @return pca plots
#' @export
plotMatPCA<-function(dataMatrix,classes,rescale=T){
  if(rescale==T){
    dataMatrix<-rescale_minus1To1(dataMatrix)
  }
  stopifnot(isMatrixValid(dataMatrix, valueRange=c(-1,1), NAsValid=FALSE))
  matpca<-stats::prcomp(dataMatrix,center=T)
  p1<-ggbiplot::ggbiplot(matpca,choices=c(1,2),labels=classes,groups=classes,
                    var.axes=F) + ggplot2::theme(legend.position = "none")
  p2<-ggbiplot::ggbiplot(matpca,choices=c(3,4),labels=classes,groups=classes,
               var.axes=F) + ggplot2::theme(legend.position = "none")
  p3<-ggbiplot::ggbiplot(matpca,choices=c(5,6),labels=classes,groups=classes,
               var.axes=F) + ggplot2::theme(legend.position = "none")
  p4<-ggbiplot::ggbiplot(matpca,choices=c(7,8),labels=classes,groups=classes,
               var.axes=F) + ggplot2::theme(legend.position = "none")
  p<-ggpubr::ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
  return(p)
}



#' Plot PCAs of matricies for different classifications
#'
#' @param k_range range of number of classes
#' @param outPath Path for out put file
#' @param outFileBase the basename of the file of the matrices
#' @return successful completion message
#' @export
plotPCAofMatrixClasses<-function(k_range, outPath, outFileBase){
  for(numClasses in k_range) {
   dataMatrix<-readRDS(paste0(outPath, "/", outFileBase, "_K",
                              numClasses, ".rds"))
   classes<-gsub("^.*__class","",rownames(dataMatrix))
   p<-plotMatPCA(dataMatrix,classes)
   p<-ggpubr::annotate_figure(p,top=ggpubr::text_grob(outFileBase,size=14))
   ggplot2::ggsave(filename=paste0(outPath,"/classPCA_",
                                   outFileBase,"_K", numClasses,".png"),
                   plot=p, device="png", width=29, height=19, units="cm")
  }
  return("PCA plotted correctly")
}



#' Plot UMAP of data matrix
#'
#' @param dataMatrix Matrix for UMAP dimentionality reduction
#' @param classes Factorised vector of classifications for all the rows
#' @param custom.settings Settings for UMAP plotting. Based on umap.defaults.
#' @return single UMAP plot
#' @export
doSingleUMAPplot<-function(dataMatrix, classes,
                           custom.settings=umap::umap.defaults){
  X1 <- X2 <- NULL
  mumap<-umap::umap(dataMatrix,config=custom.settings)
  colnames(mumap$layout)<-c("X1","X2")
  p<-ggplot2::ggplot(data=as.data.frame(mumap$layout),ggplot2::aes(x=X1,y=X2)) +
    ggplot2::geom_point(ggplot2::aes(colour=classes)) +
    ggplot2::ggtitle(paste0("UMAP with ",custom.settings$metric,
                            " distance, min_dist=", custom.settings$min_dist,
                            " spread=", custom.settings$spread))+
    ggplot2::theme(plot.title = ggplot2::element_text(size=10),
                   axis.title.x = ggplot2::element_text(size=8),
                   axis.title.y = ggplot2::element_text(size=8))
  return(p)
}




#' Plot UMAP of data matrix
#'
#' Plot Uniform Manifold Approximation and Projection (UMAP) of the reads and
#' colour by class. If matrix contains more than 2000 reads, the reads will
#' be randomly subsampled to a maximum of 2000 in order to speed up projection.
#' @param dataMatrix Matrix for UMAP dimentionality reduction
#' @param classes Vector of classifications for all the rows
#' @param rescale Rescale matrix with 0-1 values to -1 to 1 values where -1 is
#' no methylation, NAs are 0 and 1 is methylation at that position.
#' @return umap plots with different distance metrics
#' @export
plotMatUMAP<-function(dataMatrix, classes, rescale=T){
  if(rescale==T){
    dataMatrix<-rescale_minus1To1(dataMatrix)
  }
  stopifnot(isMatrixValid(dataMatrix, valueRange=c(-1,1), NAsValid=FALSE))
  dups<-duplicated(dataMatrix)
  print(paste(sum(dups)," duplicate rows removed"))
  dataMatrix<-dataMatrix[!dups,]
  readClasses<-factor(classes[!dups])

  # subsample reads if more than 2000
  if(dim(dataMatrix)[2]>2000){
    idx<-sample(1:dim(dataMatrix)[1],2000,replace=F)
    dataMatrix<-dataMatrix[idx,]
    readClasses<-readClasses[idx]
  }
  numClasses<-length(levels(readClasses))

  custom.settings = umap::umap.defaults
  #custom.settings$n_neighbors = floor(length(readClasses)/numClasses)
  #custom.settings$min_dist = 0.1 #0.1
  #custom.settings$spread=1#1

  #custom.settings$metric = "euclidean"
  #p1<-doSingleUMAPplot(dataMatrix,readClasses,custom.settings)

  custom.settings$metric = "cosine"
  p3<-doSingleUMAPplot(dataMatrix,readClasses,custom.settings)

  #p<-ggpubr::ggarrange(p3,nrow=1,ncol=1)
  return(p3)
}



#' Plot UMAP for different classifications
#'
#' @param k_range range of number of classes
#' @param outPath Path for out put file
#' @param outFileBase the basename of the file of the matrices
#' @return successful completion message
#' @export
plotUMAPofMatrixClasses<-function(k_range, outPath, outFileBase){
  for(numClasses in k_range){
    dataMatrix<-readRDS(paste0(outPath, "/", outFileBase, "_K",
                               numClasses, ".rds"))
    classes<-gsub("^.*__class","",rownames(dataMatrix))
    p<-plotMatUMAP(dataMatrix,classes)
    p<-ggpubr::annotate_figure(p,top=ggpubr::text_grob(outFileBase,size=14))
    ggplot2::ggsave(filename=paste0(outPath,"/classUMAP_",
                                    outFileBase,"_K", numClasses,".png"),
                    plot=p, device="png", width=19, height=19, units="cm")
  }
  return("UMAP plotted correctly")
}
