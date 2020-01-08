context("utils")

testthat::test_that("getMode gets mode", {
  v1<-c(2,2,3,5)
  v2<-c("b","c","a","d","a")
  v3<-c(3,3,2,2)
  testthat::expect_equal(getMode(v1), 2)
  testthat::expect_equal(getMode(v2), "a")
  testthat::expect_equal(getMode(v3), 3)
})
