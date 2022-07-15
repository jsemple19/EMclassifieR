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
#' @param figFormat format of output figures. Should be one of "png" or "pdf"
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return  allClassMeans data.frame with columns: position, methFreq, class, replicate
#' @export
runClassLikelihoodRpts<-function(dataMatrix,classes,numRepeats=20, outPath=".",
                                 xRange=c(-250,250), outFileBase="",
                                 myXlab="CpG/GpC position", featureLabel="TSS",
                                 baseFontSize=12, figFormat="png",
                                 distMetric=list(name="euclidean")){
  dataMatrix_original<-dataMatrix
  if(ncol(dataMatrix)!=ncol(classes)){
    chooseCols<-match(colnames(dataMatrix),colnames(classes))
    dataMatrix<-dataMatrix[,!is.na(chooseCols)]
    classes<-classes[,stats::na.omit(chooseCols)]
  }
  # if there are no row names
  if(is.null(rownames(classes))){
    rownames(classes)<-1:nrow(classes)
  }

  rownames(classes)<-gsub("[^0-9*]","",rownames(classes))

  numClasses<-nrow(classes)
  classVote<-matrix(nrow=nrow(dataMatrix),ncol=numRepeats)
  colnames(classVote)<-paste0("rep",1:numRepeats)
  silStats<-NULL
  allClassMeans<-NULL

  for (i in 1:numRepeats){
    lik<-classLikelihoodPerRead(dataMatrix,classes)
    classVote[,i]<-colnames(lik)[apply(lik,1,which.max)]
    classMeans<-data.frame(data.frame(dataMatrix) %>%
                    dplyr::group_split(classVote[,i], .keep=F) %>%
                    purrr::map_dfr(.f=colMeans,na.rm=T),stringsAsFactors = F)
    colnames(classMeans)<-colnames(dataMatrix)
    tmpMeans<-tidyr::gather(classMeans,key="position",value="methFreq")
    tmpMeans$class<-rep(rownames(classMeans),ncol(classMeans))
    tmpMeans$replicate<-i
    if(is.null(allClassMeans)){
      allClassMeans<-tmpMeans
    } else {
      allClassMeans<-rbind(allClassMeans,tmpMeans)
    }
  }

  classVote<-data.frame(classVote,stringsAsFactors = F)
  classVote$read<-rownames(dataMatrix)
  classVote$topClass<-NA
  classVote$topClassFreq<-NA
  classVote<-getClassVote(classVote)

  #classVote$topClassPad<-sprintf(paste0("%0",max(nchar(classVote$topClass)),"s"),trimws(classVote$topClass))
  #classVote$topClassPad<-factor(classVote$topClassPad)
  #classVote<-classVote[order(classVote$topClass),]
  classVote<-classVote[order(as.numeric(classVote$topClass)),]

  plotClassStability(classVote,outFileBase,outPath,numClasses)

  # save data with most frequent class call.
  #print("saving data with most frequent class call")
  idx<-match(row.names(dataMatrix_original),classVote$read)
  row.names(dataMatrix_original)<-paste0(rownames(dataMatrix_original),
                                        "__class",
                                  formatC(as.numeric(classVote$topClass[idx]),
                                    width=nchar(numClasses),
                                    flag=0))
  #dataOrderedByClass<-dataMatrix_original[order(rownames(classes)<-1:nrow(classes)[idx]),]
  dataOrderedByClass<-dataMatrix_original[order(
    as.numeric(classVote$topClass[idx])),]

  saveRDS(dataOrderedByClass, file=paste0(outPath, "/", outFileBase, "_K",
                                          numClasses, ".rds"))
  print("plotting final classes")
  plotFinalClasses(dataOrderedByClass, numClasses, allClassMeans,
                   outFileBase, outPath, xRange, myXlab, featureLabel,
                   baseFontSize,figFormat)

  silStats<-getSilhouetteStats(dataOrderedByClass, numClasses, outFileBase,
                               outPath, EMrep=NULL, silStats=NULL,
                               distMetric=distMetric)

  print("saving silouhette stats")
  utils::write.csv(silStats, file=paste0(outPath,"/silStats_",
                                         outFileBase,"_K",numClasses,".csv"),
                   row.names=F)

  return(dataOrderedByClass)
}




#' Get a matrix of class mean profiles from allCalssMeans table
#'
#' allClassMeans is a list of tables (one table for each number of classes). The function extracts one such table and converts it to a matrix with each row being a particular class and each column being a genomic position (in coordinates relative to an anchor point, like a TSS).
#' @param allClassMeans A list by number of classes of tables containing positions, fraction methylation, and class ID for replicate EM runs.
#' @param numClasses An integer indicating the total number of classes being used
#' @return A matrix with classID x position containing fraction methylation values for each class at each position
#' @export
getClassMeansMatrix<-function(allClassMeans,numClasses){
  position<-methFreq<-classMeans<-NULL
  dd<-allClassMeans[[numClasses]] %>%
    dplyr::group_by(position,class) %>%
    dplyr::summarize(classMeans=mean(methFreq))
  dd<-dd %>% tidyr::pivot_wider(names_from=position,
                                values_from=classMeans)
  ddMat<-as.matrix(dd[,-1])
  row.names(ddMat)<-dd$class
  ddMat<-ddMat[,order(as.numeric(colnames(ddMat)))]
  ddMat
  return(ddMat)
}




#' Pad end of matrix with missing columns
#'
#' This function compensates for shrinking of a matrix due to the application of
#' a sliding window to data. It will add the columns that are missing from
#' the edges of the matrix to the full xRange. they will take the values from
#' the first or last columns of the matrix. This might be useful for classMeans
#' derived from multi-gene windowed matrices when they have to be applied to
#' single gene matrices.
#' @param dataMat A matrix of methylation or bincount values (reads x position)
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @return A matrix with additional columns corresponding to any positions that
#' missing from the matrix relative to the xRange. These columns will take the
#' value of the firts and last columns of the original matrix
#' @export
padEndToRange<-function(dataMat, xRange=c(-250,250)){
  matStart<-min(as.numeric(colnames(dataMat)))
  matEnd<-max(as.numeric(colnames(dataMat)))
  #extra start columns
  if(matStart>xRange[1]){
    startColNames<-seq(xRange[1],matStart-1)
    startCols<-matrix(data=dataMat[,as.character(matStart)],
                      nrow=nrow(dataMat), ncol=length(startColNames))
    colnames(startCols)<-startColNames
    dataMat<-cbind(startCols,dataMat)
  }
  #extra end columns
  if(matEnd<xRange[2]){
    endColNames<-seq(matEnd+1,xRange[2])
    endCols<-matrix(data=dataMat[,as.character(matEnd)],
                    nrow=nrow(dataMat), ncol=length(endColNames))
    colnames(endCols)<-endColNames
    dataMat<-cbind(dataMat,endCols)
  }
  return(dataMat)
}

#' Collect all the class means from allClassMeans object
#'
#' AllClassMeans is a list of tables with the position in the list corresponding to the total number of classes used in clustering and the table giving the methylation frequency at different positions for each class and each replicate. This function collects all the class mean profiles for clustering into any k clusters in k_range, and returns a matrix with cluster-name-number x position
#' @param allClassMeans A list by number of classes of tables containing positions, fraction methylation, and class ID for replicate EM runs.
#' @param k_range An vector of the range of total number of classes used in each clustering
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param name Name of data set that will be added to the start of the row names.
#' @return A matrix of all the class mean profiles from the clustering into all the k clusters in k_range
#' @export
getAllClassMeansMatrix<-function(allClassMeans, k_range=4:8, xRange=c(-250,250), name="") {
  allMeansMat<-NULL
  for(k in k_range){
     classes<-getClassMeansMatrix(allClassMeans,k)
     classes<-padEndToRange(classes,xRange=c(-250,250))
     row.names(classes)<-paste0(name,"k",k,"_class",row.names(classes))
     if(is.null(allMeansMat)){
        allMeansMat<-classes
     } else {
       allMeansMat<-rbind(allMeansMat,classes)
     }
  }
  return(allMeansMat)
}



#' Combine read counts in each class from silStats files
#'
#' Function takes all region_sample listed in matTable and extracts the read
#' counts in each class from the single line silStats files that are generated
#' when clustering according to known classes.
#' @param matTable Table with list of all samples and regions which have been
#' classified according to known classes.
#' @param numClasses An integer indicating the total number of classes being used
#' @param outPath A string with the path to the directory where the output should go
#' @return data.frame of the number of reads in each class, with extra columns for
#' total number of reads, region and sample names
#' @export
combineSilStats<-function(matTable, numClasses, outPath){
  silTable<-NULL
  for(i in 1:nrow(matTable)){
    regionName=matTable$region[i]
    sampleName=matTable$sample[i]
    outFileBase=paste(sampleName, regionName, sep="_")
    silLine<-tryCatch({
      read.csv(file=paste0(outPath, "/silStats_", outFileBase,"_K",
                              numClasses, ".csv"),header=T,stringsAsFactors=F)
    },
    error=function(e){return("no silStatsFile")}
    )
    if(length(silLine)>1){
      readCounts<-silLine[grep("_reads",colnames(silLine))]
      readCounts$totalReads=sum(readCounts)
    } else {
      readCounts<-as.data.frame(rep(NA,numClasses))
      colnames(readCounts)<-paste0("class",1:numClasses,"_reads")
      readCounts$totalReads=NA
    }
    readCounts$region=regionName
    readCounts$sample=sampleName
    if(is.null(silTable)){
      silTable<-readCounts
    } else {
      silTable<-rbind(silTable,readCounts)
    }
  }
  return(silTable)
}

