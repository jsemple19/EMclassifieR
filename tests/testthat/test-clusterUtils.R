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




testthat::test_that("clusterRandomMatrices returns correct WSS", {
  dataMatrix<-matrix(rep(c(1,0,1,1,0,0.5),8),nrow=4)
  rownames(dataMatrix)<-paste0("r",1:dim(dataMatrix)[1])
  colnames(dataMatrix)<-paste0("c",1:dim(dataMatrix)[2])
  distMetric=list(name="euclidean", rescale=F)
  WSS<-clusterRandomMatrices(dataMatrix, k_range=2:3, maxB=10,
                             convergenceError=1e-6, maxIterations=100,
                             nThreads=1,setSeed=T, distMetric=distMetric)
  floor(WSS)
  truth<-data.frame(numClasses=c(2,3), meanWSS=c(10,3),
                    sumSq=c(128,16), sdWSS=c(4,1))
  testthat::expect_equal(floor(WSS),truth)
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




testthat::test_that("cosine distance object is made", {
  set.seed(200413)
  binMat<-matrix(sample(c(0,1),95,replace=T),nrow=5)
  result<-sum(cosineDist(binMat))
  testthat::expect_equal(floor(result),4)
  testthat::expect_equal(class(cosineDist(binMat)), "dist")
})



testthat::test_that("discrimination between distance metrics", {
  binMat<-matrix(rbind(c(0,0,1,1,1,1,1,0,0), c(0,0,1,0,1,0,1,0,0),
                 c(0,0,1,1,1,0,0,0,0), c(1,1,1,1,1,0,0,0,0),
                 c(0,1,1,1,1,0,0,0,0), c(0,1,1,1,1,0,0,0,1),
                 c(0,0,1,0,1,0,0,0,0), c(0,0,1,0,1,0,NA,0,0)), nrow=8)
  result1<-euclideanDist(binMat,rescale=F)
  result2<-cosineDist(binMat,valNA=0.5, rescale=F)
  result3<-cosineDist(binMat, rescale=F)
  result4<-correlationDist(binMat, valNA=0.5, rescale=T)
  testthat::expect_equal(floor(sum(result1)),39)
  testthat::expect_equal(floor(sum(result2)),7)
  testthat::expect_true(is.na(sum(result3)))
  testthat::expect_equal(floor(sum(result4)),13)
})
