context("EMmultigene")

testthat::test_that("selectReadsFromMatrix chooses correct rows", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,1,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,NA,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-selectReadsFromMatrix(dataMatrix, minReads=5, addToReadName="s1",preferBest=T)
  testthat::expect_equal(row.names(dm),c("s1__r3", "s1__r6", "s1__r4", "s1__r2",
                                         "s1__r1"))
  })


testthat::test_that("SumSqMatrixRows gets SS distance between rows of a matrix", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,NA,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,NA,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-getFullMatrix(dataMatrix,winSize=12,anchorPoint="middle")
  testthat::expect_equal(dim(dm),c(6,12))
  })



testthat::test_that("getFullMatrix expands matrix to correct dimensions", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,NA,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,NA,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-getFullMatrix(dataMatrix,winSize=12,anchorPoint="middle")
  testthat::expect_equal(dim(dm),c(6,12))
})



testthat::test_that("prepareWindows averages values correctly", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,NA,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,1,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-getFullMatrix(dataMatrix,winSize=12,anchorPoint="middle")
  dm1<-prepareWindows(dm, binSize=3, stepSize=1)
  testthat::expect_true(is.na(dm1[3,2]))
  testthat::expect_equal(dm1[1,2],0.5)
  testthat::expect_equal(dm1[2,1],0)
  testthat::expect_equal(dm1[5,9],1)
})





testthat::test_that("rescale_minus1To1 values correctly", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,NA,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,1,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-getFullMatrix(dataMatrix,winSize=12,anchorPoint="middle")
  dm1<-prepareWindows(dm, binSize=3, stepSize=1)
  dm2<-rescale_minus1To1(dm1,randomiseNAs=F)
  testthat::expect_false(is.na(dm2[3,2]))
  testthat::expect_equal(dm2[3,2],0)
  testthat::expect_equal(dm2[1,2],0)
  testthat::expect_equal(dm2[2,1],-1)
  testthat::expect_equal(dm2[5,9],1)
})



testthat::test_that("rescale_0To1 values correctly", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,NA,1,1,1,1,
                       1,1,1,0,1,1,
                       0,NA,0,1,NA,1,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-6,-3,-1,1,4,5)
  set.seed(123)
  dm<-getFullMatrix(dataMatrix,winSize=12,anchorPoint="middle")
  dm1<-prepareWindows(dm, binSize=3, stepSize=1)
  dm2<-rescale_minus1To1(dm1,randomiseNAs=F)
  dm3<-rescale_0To1(dm2)
  testthat::expect_false(is.na(dm3[3,2]))
  testthat::expect_equal(dm3[3,2],0.5)
  testthat::expect_equal(dm3[1,2],0.5)
  testthat::expect_equal(dm3[2,1],0)
  testthat::expect_equal(dm3[5,9],1)
})



testthat::test_that("countGenesPerClass works", {
  filename<-system.file("extdata", "EMres/dS02-182_multiGene_K5.rds",
              package="EMclassifieR", mustWork=TRUE)
  dataOrderedByClass<-readRDS(filename)
  counts<-countGenesPerClass(dataOrderedByClass)
  testthat::expect_equal(as.numeric(counts[[1]][1,3]),11)
  testthat::expect_true(ggplot2::is.ggplot(counts[[2]]))
})