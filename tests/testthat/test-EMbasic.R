context("EMbasic")

testthat::test_that("1 round of EM works", {
  numSamples=4
  numClasses=2
  dataMatrix<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,6),nrow=2)
  priorProb=rep(1/numClasses,numClasses)

  testthat::expect_equal(em_basic(classes,priorProb,dataMatrix)$posteriorProb,
               matrix(rep(0.5,8),nrow=4))
})


testthat::test_that("EM_basic and runEM fail if data not in 0-1 range", {
  numSamples=4
  numClasses=2
  dataMatrix<-matrix(c(0,0,0,1,1,1,2,1,0,0,0,1),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,6),nrow=2)
  priorProb=rep(1/numClasses,numClasses)
  result1<-try(em_basic(classes,priorProb,dataMatrix),silent=TRUE)
  result2<-try(runEM(dataMatrix,numClasses,1e-6,100),silent=TRUE)

  testthat::expect_equal(class(result1), "try-error")
  testthat::expect_equal(class(result2), "try-error")
})



testthat::test_that("EM_basic and runEM fail if data contains NAs", {
  numSamples=4
  numClasses=2
  dataMatrix<-matrix(c(NA,0,0,1,1,1,NA,1,0,0,0,1),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,6),nrow=2)
  priorProb=rep(1/numClasses,numClasses)
  result1<-try(em_basic(classes,priorProb,dataMatrix),silent=TRUE)
  result2<-try(runEM(dataMatrix,numClasses,1e-6,100),silent=TRUE)

  testthat::expect_equal(class(result1), "try-error")
  testthat::expect_equal(class(result2), "try-error")
})




testthat::test_that("EM_basic fails if data not in 0-1 range", {
  numSamples=4
  numClasses=2
  dataMatrix<-matrix(c(0,0,0,1,1,2,1,1,0,0,0,1),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,6),nrow=2)
  priorProb=rep(1/numClasses,numClasses)
  result<-try(em_basic(classes,priorProb,dataMatrix),silent=TRUE)
  testthat::expect_equal(class(result), "try-error")
})


testthat::test_that("1 round of EM works with fractions", {
  numSamples=5
  numClasses=3
  dataMatrix<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  classes = matrix(rep(0.5,numClasses*numSamples),nrow=numClasses)
  priorProb=rep(1/numClasses,numClasses)

  testthat::expect_equal(em_basic(classes,priorProb,dataMatrix)$posteriorProb,
               matrix(rep(1/3,numSamples*numClasses),nrow=numSamples))
})




testthat::test_that("running multiple rounds of EM  works", {
  numSamples=4
  numClasses=2
  dataMatrix<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1),nrow=numSamples,byrow=T)

  set.seed(1)
  result<-round(runEM(dataMatrix,numClasses,1e-6,100)$posteriorProb,0)
  testthat::expect_equal(result,
               matrix(c(1,0,0,1,0,1,1,0),nrow=numSamples,byrow=T))
})


testthat::test_that("running multiple rounds of EM with fractions works", {
  numSamples=5
  numClasses=3
  dataMatrix<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)

  set.seed(1)
  result<-round(runEM(dataMatrix,numClasses,1e-6,100)$posteriorProb,0)
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
  dataMatrix<-matrix(c(0,0,0,1,1,1,1,1,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-2,-1,1)
  #readClasses = apply(posteriorProb, 1, which.max)
  expectedClasses<- c("r2__class1","r3__class1","r1__class3","r4__class3","r5__class3")

  set.seed(1)
  emClass<-runEM(dataMatrix,numClasses,1e-6,100)
  orderedData<-classifyAndSortReads(dataMatrix, emClass$posteriorProb,
                                    previousClassMeans=NULL)
  results<-rownames(orderedData$data)
  testthat::expect_equal(results, expectedClasses)
})




testthat::test_that("classifying and sorting by previous class means works", {
  numClasses=3
  numSamples=6
  dataMatrix<-matrix(c(1,0.5,0,0,
                 0,0,0,0,
                 1,1,1,1,
                 1,1,1,0,
                 0,0,0,1,
                 1,0.5,0.5,0),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-paste0("r",1:numSamples)
  colnames(dataMatrix)<-c(-2,-1,1,3)
  #readClasses = apply(posteriorProb, 1, which.max)
  expectedClasses<- c("r2__class1", "r5__class1", "r3__class2", "r4__class2",
                      "r1__class3", "r6__class3")

  set.seed(1)
  emClass<-runEM(dataMatrix,numClasses,1e-6,100)
  orderedData<-classifyAndSortReads(dataMatrix, emClass$posteriorProb,
                                         previousClassMeans=NULL)
  results<-rownames(orderedData$data)
  previousClassMeans<-orderedData$classMeans
  orderedData1<-classifyAndSortReads(dataMatrix, emClass$posteriorProb,
                                    previousClassMeans=previousClassMeans)
  results1<-rownames(orderedData1$data)
  #classMeans<-orderedData1$classMeans
  testthat::expect_equal(results, expectedClasses)
  testthat::expect_equal(results, results1)
})




testthat::test_that("plotting from single gene returns a plot" , {
  numSamples<-5
  dataMatrix<-matrix(c(1,1,1,1,1,0,0,0,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-c("r2__class1","r3__class1",
                    "r1__class3","r4__class3","r5__class3")
  colnames(dataMatrix)<-c(-200,-100, 200)

  #testthat::expect_null(plotClassesSingleGene(matrix(c(0)),performSilhouette=F),)
  p<-plotClassesSingleGene(dataMatrix, xRange=c(-250,250),
                        title="Reads by classes",
                        myXlab="CpG/GpC position",
                        featureLabel="TSS",
                        baseFontSize=12)
  testthat::expect_equal(class(p),c("gg","ggplot"))
})



testthat::test_that("silhouette plot works" , {
  numSamples<-5
  dataMatrix<-matrix(c(1,1,1,1,1,0,0,0,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-c("r2__class1","r3__class1","r1__class3","r4__class3","r5__class3")
  colnames(dataMatrix)<-c(-200,-100, 200)

  sp<-silhouettePlot(dataOrderedByClass=dataMatrix, numClasses=3, outFileBase="blah")
  testthat::expect_equal(class(sp$stats),"data.frame")
  testthat::expect_equal(class(sp$plotObject),"silhouette")
})





testthat::test_that("plotting class means returns a plot" , {
  numSamples<-5
  numClasses<-3
  dataMatrix<-matrix(c(1,1,1,1,1,0,0,0,0,0,0,1,0.5,0.5,0.5),nrow=numSamples,byrow=T)
  rownames(dataMatrix)<-c("r2__class1","r3__class1",
                    "r1__class3","r4__class3","r5__class3")
  colnames(dataMatrix)<-c(-200,-100, 200)
  emClass<-runEM(dataMatrix,numClasses,1e-6,100)

  p<-plotClassMeans(emClass$classes, xRange=c(-250,250), facet=TRUE, title="Class means",
                    myXlab="CpG/GpC position",featureLabel="TSS", baseFontSize=12)
  testthat::expect_equal(class(p),c("gg","ggplot"))
})




testthat::test_that("plotting smoothed class means returns a plot" , {
  df<-data.frame(position=c(-2,-2,-2,-1,-1,-1,2,2,2),
                 methFreq=rep(c(0.2,0.5,0.8),3)+runif(n=9,min=0.05,max=0.15),
                 class=factor(rep(1:3,3)),
                 replicate=factor(rep(1:3,each=3)),stringsAsFactors=FALSE)

  p<-plotSmoothedClassMeans(df, xRange=c(-2,2))
  testthat::expect_equal(class(p),c("gg","ggplot"))
})




testthat::test_that("runEMrepeats class mean output is correct" , {
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

  set.seed(1)
  df<-runEMrepeats(dataMatrix, numClasses=numClasses, EMrepeats=3, xRange=c(-3,3))
  testthat::expect_equal(dim(df),c(36,4))
  testthat::expect_equal(df[3,3],"1")
})





testthat::test_that("getClassVote returns table with most frequent class", {
  df<-data.frame(read=c("r1","r2","r3"),rep1=c(3,2,2),rep2=c(3,1,2),rep3=c(3,2,2),rep4=c(3,2,2))
  df1<-getClassVote(df)
  testthat::expect_equal(df1$topClass,c(3,2,2))
  testthat::expect_equal(df1$topClassFreq,c(1.00,0.75,1.00))
})



