## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressMessages(library(EMclassifieR))

## -----------------------------------------------------------------------------
tablePath="csv/MatrixLog_relCoord_ampTSS.csv"
matTable<-read.csv(system.file("extdata", tablePath, package = "EMclassifieR",
                               mustWork=TRUE), stringsAsFactors=F)
head(matTable)

## -----------------------------------------------------------------------------
i=1
dataMatrix<-readRDS(system.file("extdata", matTable$filename[i], package = "EMclassifieR", mustWork=TRUE))
dim(dataMatrix)
dataMatrix<-removeAllNArows(dataMatrix)
dim(dataMatrix)

## ----eval=T-------------------------------------------------------------------

k_range = 2:8      # Number of classes to be found
maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
convergenceError = 10e-6
numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
xRange=c(-250,250)
maxB=50 # Number of randomised matrices to generate
outPath="./EMres"
maxTasks=4
taskId=1
nThreads=4
setSeed=FALSE

#split table indecies into nTasks number of groups
taskSubList<-split(1:nrow(matTable),sort(1:nrow(matTable)%%maxTasks))

set.seed(200413)
for (i in taskSubList[[taskId]]){
#for (i in 1:nrow(matTable)) {
  regionName=matTable$region[i]
  sampleName=matTable$sample[i]
  outFileBase=paste(sampleName, regionName, sep="_")
  print(paste("Clustering", outFileBase))
  dataMatrix<-readRDS(system.file("extdata", matTable$filename[i], 
                              package = "EMclassifieR", mustWork=TRUE))
  dim(dataMatrix)
  dataMatrix<-removeAllNArows(dataMatrix)
  dim(dataMatrix)
  
  allClassMeans<-tryCatch( 
    {
      print("running EM for a range of class numbers")
      runEMrangeClassNum(dataMatrix, k_range, convergenceError, maxIterations,
                     EMrepeats=numRepeats, outPath=outPath, xRange=xRange, 
                     outFileBase=paste(sampleName, regionName, sep="_"),
                     doIndividualPlots=FALSE)
    },
    error=function(e){"Matrix not valid"}
    )
      
  if(is.list(allClassMeans)){
     	saveRDS(allClassMeans,paste0(outPath,"/allClassMeans_",outFileBase,".rds"))
  } else {
     	print(allClassMeans) # error message
  }

  clustMetrics<-tryCatch( {
	print("plotting clustering metrics for a range of class sizes")
	plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
		maxIterations, outPath, outFileBase, EMrep=NULL, nThreads=nThreads, 
		setSeed=setSeed)
  },
  error=function(e){"Matrix not valid"}
  )
  if(length(clustMetrics)==1) {
    print(clustMetrics)
  }
}


