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

## ----eval=F-------------------------------------------------------------------
#  
#  k_range = 2:8      # Number of classes to be found
#  maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
#  convergenceError = 10e-6
#  numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
#  xRange=c(-250,250)
#  maxB=100 # Number of randomised matrices to generate
#  outPath="./EMres"
#  
#  #nThreads<-detectCores()-1
#  
#  set.seed(20200413)
#  for (i in 1:nrow(matTable)){
#    regionName=matTable$region[i]
#    sampleName=matTable$sample[i]
#    outFileBase=paste(sampleName, regionName, sep="_")
#    print(paste("Clustering", outFileBase))
#    dataMatrix<-readRDS(system.file("extdata", matTable$filename[i],
#                                package = "EMclassifieR", mustWork=TRUE))
#    dim(dataMatrix)
#    dataMatrix<-removeAllNArows(dataMatrix)
#    dim(dataMatrix)
#  
#    tryCatch(
#      {
#        runEMrangeClassNum(dataMatrix, k_range, convergenceError, maxIterations,
#                       repeats=numRepeats, outPath=outPath, xRange=xRange,
#                       outFileBase=paste(sampleName, regionName, sep="_"),
#                       doIndividualPlots=FALSE)
#  
#        plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
#                         maxIterations, outPath, outFileBase)
#      },
#      error=function(e){"Matrix not valid"}
#    )
#  }
#  

