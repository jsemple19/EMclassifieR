context("EMbasic")

testthat::test_that("1 round of EM works", {
  numSamples=4
  numClasses=2
  data<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,6),nrow=2)
  priorProb=rep(1/numClasses,numClasses)

  testthat::expect_equal(em_basic(classes,priorProb,data)$posteriorProb,
               matrix(rep(0.5,8),nrow=4))
})


testthat::test_that("1 round of EM works with fractions", {
  numSamples=5
  numClasses=3
  data<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,numClasses*numSamples),nrow=numClasses)
  priorProb=rep(1/numClasses,numClasses)

  testthat::expect_equal(em_basic(classes,priorProb,data)$posteriorProb,
               matrix(rep(1/3,numSamples*numClasses),nrow=numSamples))
})




testthat::test_that("running multiple rounds of EM  works", {
  numSamples=4
  numClasses=2
  data<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1),nrow=numSamples,byrow=T)

  set.seed(1)
  result<-round(runEM(data,numClasses,1e-6,100)$posteriorProb,0)
  testthat::expect_equal(result,
               matrix(c(1,0,0,1,0,1,1,0),nrow=numSamples,byrow=T))
})


testthat::test_that("running multiple rounds of EM with fractions works", {
  numSamples=5
  numClasses=3
  data<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)

  set.seed(1)
  result<-round(runEM(data,numClasses,1e-6,100)$posteriorProb,0)
  testthat::expect_equal(result,
               matrix(c(0,0,1,1,0,0,1,0,0,0,0,1,0,0,1),nrow=numSamples,byrow=T))
})




testthat::test_that("ordering by previous class means works", {
  numClasses=2
  classMeans<-matrix(c(1,1,0.5,0,0,0.5),nrow=numClasses,byrow=T)
  prev_classMeans<-classMeans[c(2,1),]

  testthat::expect_equal(order_by_prev_cluster(numClasses,prev_classMeans,classMeans),
               c(2,1))
})



testthat::test_that("classifying and sorting reads works", {
  numClasses=3
  numSamples=5
  data<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  row.names(data)<-paste0("r",1:numSamples)
  #readsClasses = apply(posteriorProb, 1, which.max)
  expectedClasses<- c("r2__class1","r3__class1","r1__class3","r4__class3","r5__class3")

  set.seed(1)
  posteriorProb<-runEM(data,numClasses,1e-6,100)$posteriorProb
  results<-row.names(classifyAndSortReads(data, posteriorProb, previousClassMeans=NULL))
  testthat::expect_equal(results, expectedClasses)
})



testthat::test_that("plotting from single gene works" , {
  numSamples<-5
  data<-matrix(c(1,1,1,1,1,0,0,0,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  rownames(data)<-c("r2__class1","r3__class1","r1__class3","r4__class3","r5__class3")
  colnames(data)<-c(-200,-100, 200)

  #testthat::expect_null(plotClassesSingleGene(matrix(c(0)),performSilhouette=F),)
  testthat::expect_is(plotClassesSingleGene(data, performSilhouette=T,
                                            xlim=c(-250,250),
                                            title="Reads by classes",
                                            myXlab="CpG/GpC position",
                                            featureLabel="TSS",
                                            baseFontSize=12),"ggplot")
})
