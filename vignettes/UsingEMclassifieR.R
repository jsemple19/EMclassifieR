## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressMessages(library(EMclassifieR))

## -----------------------------------------------------------------------------
tablePath="csv/MatrixLog_relCoord_ampTSS.csv"
matTable<-read.csv(system.file("extdata", tablePath, package = "EMclassifieR",
                               mustWork=TRUE), stringsAsFactors=F)
head(matTable)

