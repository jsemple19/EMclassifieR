context("clusterUtils")

testthat::test_that("SumSqMatrixRows gets SS distance between rows of a matrix", {
  dataMatrix<-matrix(c(1,2,3,2,3,4,1,4,2,3,1,1),nrow=4)
  testthat::expect_equal(SumSqMatrixRows(dataMatrix),43)
})



testthat::test_that("withinClusterSS gets total within cluster SS distance", {
  dataMatrix<-matrix(c(1,2,3,2,3,4,1,4,2,3,1,1),nrow=4)
  classes<-c(2,1,1,2)
  testthat::expect_equal(withinClusterSS(dataMatrix,classes),17)
})


testthat::test_that("randomiseMatrixRows returns a matrix randomised by row", {
  dataMatrix<-matrix(c(1,2,3,2,3,4,1,4,2,3,1,1),nrow=4,byrow=T)
  postRand<-matrix(c(1,2,3,2,3,4,2,1,4,1,3,1),nrow=4,byrow=T)
  set.seed(1)
  testthat::expect_equal(randomiseMatrixRows(dataMatrix),postRand)
})




testthat::test_that("clusterRandomiseMatrices returns correct WSS", {
  numSamples=6
  numClasses=2
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                 0,0,0,0,0,1,
                 1,1,1,1,1,1,
                 1,1,1,0,1,1,
                 0,0,0,1,0,0,
                 1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-3,-2,-1,1,2,3)
  randWSSsd<-c(3.0,1.2)
  set.seed(1)
  WSS<-clusterRandomMatrices(dataMatrix, k_range=2:3, maxB=20,
                             convergenceError=1e-6, maxIterations=100)
  testthat::expect_equal(round(WSS$sdWSS,1),randWSSsd)
})


testthat::test_that("isMatrixValid returns TRUE for valid matrix", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1,
                       1,1,1,1,1,1,
                       1,1,1,0,1,1,
                       0,0,0,1,0,0,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-3,-2,-1,1,2,3)
  testthat::expect_equal(isMatrixValid(dataMatrix,valueRange=c(0,1)),TRUE)
  testthat::expect_equal(isMatrixValid(dataMatrix,valueRange=c(1,2)),FALSE)
})



testthat::test_that("isMatrixValid returns FALSE for matrix with NAs", {
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,NA,0,0,
                       0,0,0,0,0,1,
                       1,1,1,1,1,1,
                       1,1,1,0,1,1,
                       0,0,0,1,0,0,
                       1,0.5,1,0,1,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-3,-2,-1,1,2,3)
  NAsBad<-isMatrixValid(dataMatrix,valueRange=c(0,1),NAsValid=FALSE)
  testthat::expect_equal(NAsBad,FALSE)
  testthat::expect_equal(isMatrixValid(dataMatrix,valueRange=c(0,1),NAsValid=TRUE),TRUE)
})



testthat::test_that("isMatrixValid returns FALSE for empty or non-matrices", {
  numSamples=2
  dataMatrix<-matrix(c(1,0.5,0,0,0,0,
                       0,0,0,0,0,1),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-3,-2,-1,1,2,3)

  badMat<-isMatrixValid(dataMatrix[-c(1:2),],valueRange=c(0,1),NAsValid=FALSE)
  notMat<-isMatrixValid(c(0,0.2,1),valueRange=c(0,1),NAsValid=FALSE)
  testthat::expect_equal(badMat,FALSE)
  testthat::expect_equal(notMat,FALSE)
})
