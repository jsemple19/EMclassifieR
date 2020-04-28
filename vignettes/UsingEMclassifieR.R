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

