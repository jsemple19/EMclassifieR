context("preProcessMethMats")

testthat::test_that("NA rows are removed", {
  mat<-matrix(c(1,2,3,NA,4,5),nrow=3,byrow=T)

  testthat::expect_equal(removeAllNArows(mat),matrix(c(1,2,4,5),nrow=2,byrow=T))
})



testthat::test_that("NAs are added for missing columns", {
  mat<-matrix(c(1:6),nrow=2,byrow=T)
  colnames(mat)<-c(-2,1,2)
  rownames(mat)<-c("a","b")
  endMat<-matrix(c(1,NA,2,3,4,NA,5,6),nrow=2,byrow=T)
  colnames(endMat)<-c(-2,-1,1,2)
  rownames(endMat)<-c("a","b")

  testthat::expect_equal(padMissingColWithNAs(mat,colRange=c(-2,2)),
                                              endMat)
})



testthat::test_that("0 is removed from the vector", {
  vec<-c(-2:2)

  testthat::expect_equal(remove0(vec),c(-2,-1,1,2))
})


test_that("does not fail when 0 is absent", {
  vec<-c(-2,-1,1,2)

  expect_equal(remove0(vec),c(-2,-1,1,2))
})