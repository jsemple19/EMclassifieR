# classifyByKnownClasses.R
# Functions for classifying reads by known class profiles

#' Get likelihood of different classes for each read
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param classes A matrix of classes x position giving the mean frequency of methylation at each position for each class
#' @return likelihood of different classes for each read
#' @export
classLikelihoodPerRead<-function(dataMatrix,classes){
  dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
  stopifnot(isMatrixValid(dataMatrix, valueRange=c(0,1),NAsValid=FALSE))
  fractionIdx<-dataMatrix > 0 & dataMatrix < 1
  binaryData<-dataMatrix
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx), size=1,
                                           prob=dataMatrix[fractionIdx])
  }

  numClasses=dim(classes)[1]         # Number of classes
  numSamples=dim(binaryData)[1]            # Number of samples
  l=matrix(nrow=numSamples, ncol=numClasses)  # log likelihood (theta|data)

  # calculate probability of seeing this pattern given the classes (i.e. likelihood)
  for(i in 1:numSamples) {
    for (j in 1:numClasses) {
      # Bernoulli (size = 1): simply multiplying by class probability
      l[i,j] = sum(stats::dbinom(x = binaryData[i,], size = 1, prob = classes[j,],
                                 log=T))
    }
  }

  # normalise and convert to probabilities per sample
  for(i in 1:numSamples) {
    l[i,]<-exp(l[i,])/sum(exp(l[i,]))
  }
  colnames(l)<-rownames(classes)
  rownames(l)<-rownames(dataMatrix)
  return(l)
}





#' Run several repeats classifying reads into classes
#'
#' Reads are classified into different classes based on likelihood estimation using the class means. Due to the need for randomisation when binarisation of the methylation matrix, the classification is carried out repeatedly and then the most frequently assigned class is chosen.
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param classes A matrix of classes x position giving the mean frequency of methylation at each position for each class
#' @param numRepeats An integer indicating the number of times to repeat the clustering (default=10)
#' @param outPath A string with the path to the directory where the output should go
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced (default is "")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return  allClassMeans data.frame with columns: position, methFreq, class, replicate
#' @export
runClassLikelihoodRpts<-function(dataMatrix,classes,numRepeats=20, outPath=".",
                                 xRange=c(-250,250), outFileBase="",
                                 myXlab="CpG/GpC position", featureLabel="TSS",
                                 baseFontSize=12, distMetric=list(name="euclidean")){
  dataMatrix_original<-dataMatrix
  if(ncol(dataMatrix)!=ncol(classes)){
    chooseCols<-match(colnames(dataMatrix),colnames(classes))
    dataMatrix<-dataMatrix[,!is.na(chooseCols)]
    classes<-classes[,stats::na.omit(chooseCols)]
  }

  numClasses<-max(as.numeric(rownames(classes)))
  classVote<-matrix(nrow=nrow(dataMatrix),ncol=numRepeats)
  colnames(classVote)<-paste0("rep",1:numRepeats)
  silStats<-NULL
  allClassMeans<-NULL

  for (i in 1:numRepeats){
    lik<-classLikelihoodPerRead(dataMatrix,classes)
    classVote[,i]<-apply(lik,1,which.max)
    classMeans = stats::aggregate(dataMatrix, by=list(classVote[,i]), FUN=mean)[-1]
    tmpMeans<-tidyr::gather(classMeans,key="position",value="methFreq")
    tmpMeans$class<-rep(rownames(classMeans),ncol(classMeans))
    tmpMeans$replicate<-i
    if(is.null(allClassMeans)){
      allClassMeans<-tmpMeans
    } else {
      allClassMeans<-rbind(allClassMeans,tmpMeans)
    }
  }

  classVote<-data.frame(classVote)
  classVote$read<-rownames(dataMatrix)
  classVote$topClass<-NA
  classVote$topClassFreq<-NA
  classVote<-getClassVote(classVote)

  classVote$topClass<-factor(classVote$topClass)
  plotClassStability(classVote,outFileBase,outPath,numClasses)

  # save data with most frequent class call.
  #print("saving data with most frequent class call")
  idx<-match(row.names(dataMatrix_original),classVote$read)
  row.names(dataMatrix_original)<-paste0(rownames(dataMatrix_original),"__class",
                                classVote$topClass[idx])
  dataOrderedByClass<-dataMatrix_original[order(classVote$topClass[idx]),]


  saveRDS(dataOrderedByClass, file=paste0(outPath, "/", outFileBase, "_K",
                                          numClasses, ".rds"))
  print("plotting final classes")
  plotFinalClasses(dataOrderedByClass, numClasses, allClassMeans,
                   outFileBase, outPath, xRange, myXlab, featureLabel,
                   baseFontSize)

  silStats<-getSilhouetteStats(dataOrderedByClass, numClasses, outFileBase,
                               outPath, EMrep=NULL, silStats=NULL,
                               distMetric=distMetric)

  print("saving silouhette stats")
  utils::write.csv(silStats, file=paste0(outPath,"/silStats_",
                                         outFileBase,"_K",numClasses,".csv"),
                   row.names=F)

  return(dataOrderedByClass)
}
