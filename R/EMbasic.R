#' Expectation Maximisation algorithm
#'
#' @param classes A matrix containing the classes to be optimised. c[i,j] is the expected methylation value or bincount value of class i at position j.
#' @param priorProb A vector defining the prior probabilities of each class.
#' @param dataMatrix  A matrix containing the sample. dataMatrix [i,j] is the methylation value or bincount value of sample i at position j.
#' @return A list of three items:
#'             1) classes: a matrix with classes being optimised (class x position)
#'             2) priorProb: a vector with the prior probabilities of each class
#'             3) posteriorProb: a matrix of probabilites of each sample belonging
#'             to a particular class (samples x class)
#' @export
em_basic <- function(classes,priorProb,dataMatrix) {
  # To deal with fractions in methlyation dataMatrix (when opposite strand do not agree on
  # methylation status), the methylation status will be applied at random with
  # probability given by the fraction.
  dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
  fractionIdx<-dataMatrix > 0 & dataMatrix < 1
  stopifnot(isMatrixValid(dataMatrix,NAsValid=FALSE))
  binaryData<-dataMatrix
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx),size=1,prob=dataMatrix[fractionIdx])
  }

  numClasses=dim(classes)[1]         # Number of classes
  numSamples=dim(binaryData)[1]            # Number of samples
  l=matrix(nrow=numSamples, ncol=numClasses)  # log likelihood (theta|data)
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability by which each class occurs

  # E step
  for(i in 1:numSamples) {
    for (j in 1:numClasses) {
      # Bernoulli (size = 1): simply multiplying by class probability
      l[i,j] = sum(stats::dbinom(x = binaryData[i,], size = 1, prob = classes[j,], log = T))
      # not sure about implementing the beta distribution
      #l[i,j] = sum(dbeta(x = dataMatrix[i,], shape1 = classes[j,], shape2 = classes[j,], log = T))
    }
  }

  # M step
  for(i in 1:numSamples) {
    posteriorProb[i,] = priorProb*exp(l[i,]-max(l[i,]))
    posteriorProb[i,] = posteriorProb[i,]/sum(posteriorProb[i,])
  }
  priorProb = colMeans(posteriorProb)
  # The expected bincounts for each class are computed at once by means of a matrix multiplication.
  classes = (t(posteriorProb) %*% binaryData)/colSums(posteriorProb)

  # classes, priorProb, posteriorProb are exported for use outside this function.
  # The output variable posterioProb is a matrix containing the samples’ posterior
  # probabilities of belonging to particular classes (rows correspond to samples,
  # columns to classes).
  #classes <<- classes; priorProb <<- priorProb; posteriorProb <<- posteriorProb;
  return(list(classes=classes,priorProb=priorProb,posteriorProb=posteriorProb))
}



#' Expectation Maximisation algorithm
#'
#' @param classes A matrix containing the classes to be optimised. c[i,j] is the expected methylation value or bincount value of class i at position j.
#' @param priorProb A vector defining the prior probabilities of each class.
#' @param dataMatrix  A matrix containing the sample. dataMatrix [i,j] is the methylation value or bincount value of sample i at position j.
#' @return A list of three items:
#'             1) classes: a matrix with classes being optimised (class x position)
#'             2) priorProb: a vector with the prior probabilities of each class
#'             3) posteriorProb: a matrix of probabilites of each sample belonging
#'             to a particular class (samples x class)
#' @export
em_classify <- function(classes,priorProb,dataMatrix) {
  # To deal with fractions in methlyation dataMatrix (when opposite strand do not agree on
  # methylation status), the methylation status will be applied at random with
  # probability given by the fraction.
  dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
  fractionIdx<-dataMatrix > 0 & dataMatrix < 1
  stopifnot(isMatrixValid(dataMatrix,NAsValid=FALSE))
  binaryData<-dataMatrix
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx),size=1,prob=dataMatrix[fractionIdx])
  }

  numClasses=dim(classes)[1]         # Number of classes
  numSamples=dim(binaryData)[1]            # Number of samples
  l=matrix(nrow=numSamples, ncol=numClasses)  # log likelihood (theta|data)
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability by which each class occurs

  # E step
  for(i in 1:numSamples) {
    for (j in 1:numClasses) {
      # Bernoulli (size = 1): simply multiplying by class probability
      l[i,j] = sum(stats::dbinom(x = binaryData[i,], size = 1, prob = classes[j,], log = T))
      # not sure about implementing the beta distribution
      #l[i,j] = sum(dbeta(x = dataMatrix[i,], shape1 = classes[j,], shape2 = classes[j,], log = T))
    }
  }
}








############################### EM HELPERS #########################################


#' Check if p value has converged
#'
#' @param p_diff The difference in likelihood betwween this round of optimisation and the previous round.
#' @param p_diff_prev The difference in likelihood between the previous round of optimisation and the one before that.
#' @param error The maximum tolerated error
#' @return A TRUE or FALSE value indicating if the difference between successive rounds of optimisation is below the error threshold (i.e. converged)
#' @export
p_converged <- function(p_diff, p_diff_prev, error) {
  return(abs(p_diff-p_diff_prev) < error)
}


#' Get absolute value of difference between two p values
#'
#' @param p Current probability
#' @param p_prev Probability from last round of optimisation
#' @return The absolute value of the change in probability
#' @export
p_difference <- function(p, p_prev) {
  return(sum(abs(p-p_prev)))
}



#' Run EM iteratively to convergence
#'
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param printProgress Print messages showing progress of convergence
#' @return  list of three items:
#'             1) classes: a matrix with optimised classes (class x position)
#'             2) priorProb: a vector with the prior probabilities of each class
#'             3) posteriorProb: a matrix of probabilites of each sample belonging
#'             to a particular class (samples x class)
#' @export
runEM<-function(dataMatrix, numClasses, convergenceError=1e-6, maxIterations=100,
                printProgress=FALSE) {
  ################################# RUN EM #########################################
  #
  # Complete algorithm for partitioning with random seeds

  # To deal with fractions in methlyation dataMatrix (when opposite strands do
  # not agree on methylation status), the methylation status will be applied at
  # random with probability given by the fraction. NAs will be recoded to 0.5
  # (equal probability of 0 or 1 at that position)
  fractionIdx<-dataMatrix > 0 & dataMatrix < 1
  dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
  stopifnot(isMatrixValid(dataMatrix, NAsValid=FALSE))
  binaryData<-dataMatrix
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx), size=1,
                                           prob=dataMatrix[fractionIdx])
  }

  numSamples=dim(binaryData)[1]            # Number of samples
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability each sample (read) belongs to a particular class

  for(i in 1:numClasses) {
    # Samples are randomly assigned probabilities (weights) for each class.
    posteriorProb[,i] = stats::rbeta(numSamples,numSamples**-0.5,1)
  }
  # these probabilities are used to generate expected bincount vectors for each class.
  # The probabilistic class assignment makes sure that classes will be free of zero values.
  # (Initial zero values are undesirable as they cannot be changed during EM.).
  classes = (t(posteriorProb) %*% binaryData)/colSums(posteriorProb)

  # a vector defining the prior probabilities of each class.
  priorProb=rep(1/numClasses,numClasses)


  e = convergenceError * length(posteriorProb) # scale convergence error to size of p

  # first iteration
  classes[classes>1]=1 # dbern error when probability > 1
  p_prev = posteriorProb
  results<-em_basic(classes,priorProb,binaryData)
  classes=results$classes
  posteriorProb=results$posteriorProb
  priorProb=results$priorProb
  p_diff = p_difference(posteriorProb, p_prev)
  p_diff_prev = 0

  i = 2
  if (printProgress) {
    print(p_converged(p_diff, p_diff_prev, e) & i < maxIterations)
  }
  while (!p_converged(p_diff, p_diff_prev, e) & i < maxIterations) {
    if (printProgress) {
      print(paste("i:", i, " , sum:", p_diff, " , step diff:", abs(p_diff-p_diff_prev)))
    }
    classes[classes>1]=1 # dbern error when probability > 1
    p_prev = posteriorProb
    p_diff_prev = p_diff

    results<-em_basic(classes,priorProb,binaryData)
    classes=results$classes
    posteriorProb=results$posteriorProb
    priorProb=results$priorProb
    p_diff = p_difference(posteriorProb, p_prev)

    i = i+1
  }
  #if (printProgress) {
  print(paste("converged:", p_converged(p_diff, p_diff_prev, e) & i < maxIterations,
              "iterations: ",i))
  #}
  return(list(classes=classes,priorProb=priorProb,posteriorProb=posteriorProb))
}



#' Order classes
#'
#' Order classes either according to the order of previous classes or order by similarity using hclust
#' @param numClasses An integer with the number of classes learned
#' @param classMeans The mean profile of all reads in their classes
#' @param prev_classMeans The mean profile of all reads in their classes from a previous run
#' @return Returns a vector with the order of the classes
#' @export
order_by_prev_cluster <- function(numClasses, classMeans, prev_classMeans) {
  cm<-classMeans
  pr = prev_classMeans
  ord = numeric(length = numClasses)
  for (i in 1:numClasses) {
    # compare cluster means by minimum sum of squares
    pos = which.min(apply(classMeans, 1, function(row) sum((row-pr[i,])^2) ))
    ord[i] = pos
    cm[pos,] = Inf
  }
  return(ord)
}



#' Classify and sort reads into their classes
#'
#' Classify reads by their posterior probability of belonging to a specific class
#' Then sort the classes by using similarity to the mean profile of previous classes.
#' If the means of previous classes was not provided, hclust is used to cluster classes by their similarity.
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param posteriorProb posteriorProb: a matrix of probabilites of each sample belonging to a particular class (samples x class)
#' @param previousClassMeans A matrix of the class means from a previous round of clustering
#' @return Returns a matrix with the reads classified (__classX is appended to the read name), and the classes are sorted.
#' @export
classifyAndSortReads<-function(dataMatrix,posteriorProb,previousClassMeans=NULL) {
  ###################### POST-EM DATA SORTING ###########################################
  #
  # assign classes to reads according to the highest class probability
  numClasses=ncol(posteriorProb)
  readClasses = apply(posteriorProb, 1, which.max)
  #readsTable = table(readClasses)
  classMeans = stats::aggregate(dataMatrix, by=list(readClasses), FUN=mean)[-1]

  if (!is.null(previousClassMeans)) {
    #print("orderByPreviousClusters")
    # order the classes by comparing the class means to the previous clustering
    classOrder = order_by_prev_cluster(numClasses, classMeans, previousClassMeans)
    classMeans = classMeans[classOrder, ]
    rownames(classMeans)<-c(1:numClasses)[c(1:numClasses) %in%
                                            unique(readClasses)]
  } else {
    #print("orderByHClust")
    # order the classes by similarity (class means)
    classOrder = stats::hclust(stats::dist(classMeans))$order
    classMeans = classMeans[classOrder, ]
    rownames(classMeans)<-c(1:numClasses)[c(1:numClasses) %in%
                                            unique(readClasses)]
  }
  #print(classOrder)
  # change name of classes to match the order
  readClasses<-factor(readClasses)
  # remap the values in the vector
  rc<-plyr::mapvalues(readClasses,
                      from = levels(readClasses)[classOrder],
                      to = levels(readClasses))
  # change the level order
  rc1<-factor(rc,levels=c(1:numClasses))
  rownames(dataMatrix)<-paste(rownames(dataMatrix),rc1,sep="__class")
  ord = order(match(rc1,1:numClasses))
  dataOrderedByClass = dataMatrix[ord,]

  return(list(dataMatrix=dataOrderedByClass, classMeans=classMeans))
}



#' Plot reads by class for a single gene
#'
#' Create a single molecule plot of the reads sorted by class.
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param title A title for the plot (default is "Reads by classes")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @param segmentSize Length of colour segment denoting methylation site
#' @param colourChoice A list of colours for colour pallette. Must include
#' values for "low", "mid", "high" and "bg" (background) and "lines".
#' @return Returns a ggplot2 object of a single molecule plot sorted by classes
#' @export
plotClassesSingleMolecule<-function(dataOrderedByClass,
                                xRange=c(-250,250), title="Reads by classes",
                                myXlab="CpG/GpC position",
                                featureLabel="TSS", baseFontSize=12,
                                segmentSize=5,
                                colourChoice=list(low="blue", mid="white",
                                                  high="red", bg="white",
                                                  lines="grey80")){
  position <- methylation <- molecules <- readNumber <- Class <-NULL
  readClasses <- paste0("class",sapply(strsplit(rownames(dataOrderedByClass),split="__class"),"[[",2))
  classOrder <- unique(readClasses)
  readNames<-sapply(strsplit(rownames(dataOrderedByClass),split="__class"),"[[",1)
  readsTable <- table(readClasses)
  #classMeans = stats::aggregate(dataMatrix, by = list(readClasses), FUN = mean)[-1]
  # the horizontal red lines on the plot
  classBorders <- utils::head(cumsum(readsTable[classOrder]), -1)+0.5
  df<-as.data.frame(dataOrderedByClass,stringsAsFactors=F)
  df1<-data.frame(read=readNames, readNumber=1:length(readNames),
                  Class=factor(readClasses), stringsAsFactors=F)

  #######################################################################################
  reads<-row.names(df)
  d<-tidyr::gather(df,key=position,value=methylation)
  d$molecules<-seq_along(reads)
  #d$methylation<-as.character(d$methylation)
  if(grepl("_",d$position[1])){
    starts<-as.numeric(sapply(strsplit(d$position,"_"),"[[",1))
    ends<-as.numeric(sapply(strsplit(d$position,"_"),"[[",2))
    d$position<-floor(rowMeans(cbind(starts,ends)))
  } else{
    d$position<-as.numeric(d$position)
  }
  p<-ggplot2::ggplot(d,ggplot2::aes(x=position,y=molecules)) +
    ggplot2::geom_tile(ggplot2::aes(width=segmentSize,fill=methylation),alpha=0.8) +
    ggplot2::scale_fill_gradient2(low=colourChoice$low, mid=colourChoice$mid,
                                  high=colourChoice$high,
                                  midpoint=0.5, na.value="transparent",
                                 breaks=c(0,1), labels=c("protected","accessible"),
                                 limits=c(0,1), name="dSMF\n\n") +
    #ggplot2::scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white", labels=c("protected","accessible"),name="dSMF") +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background=ggplot2::element_rect(fill=colourChoice$bg),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position="right", legend.box = "vertical",
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width=ggplot2::unit(0.3,"cm")) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("Single molecules") +
    ggplot2::xlim(xRange[1],as.numeric(xRange[2])+10)
    # add line for TSS
  p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=0,
                                              ymax=length(reads) +
                                                max(3, 0.04*length(reads))),
                                              color="grey80")+
      ggplot2::annotate(geom="text", x=1,
                        y=-max(2,0.03*length(reads)),
                        label=featureLabel,color="grey20")
  # add lines separating classes
  p<-p+ggplot2::geom_hline(yintercept=classBorders,colour=colourChoice$lines)
  # add color bar for classes
  p<-p+ggplot2::geom_segment(data=df1, mapping=ggplot2::aes(x=(xRange[2]+10),
                                                             y=readNumber-0.5,
                                                             xend=(xRange[2]+10),
                                                             yend=readNumber+0.5,
                                                             colour=Class),
                             size=5) +
      ggplot2::geom_vline(xintercept=xRange[2],colour="grey80")
  return(p)
}



#' Plot silhouette plot to evaluate classification
#'
#' Create a silhouette plot to evaluate the quality of classification, and return some basic parameters about the classification
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param numClasses An integer with the number of classes learned
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced
#' @param EMrep An integer indicating which EM repeat this is
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return A list of two items: a dataframe with silhouette stats and a silhouette plot object. The dataframe contains mean and SD of silhouette width overall, and per class, as well as number of reads per class
#' @export
silhouettePlot<-function(dataOrderedByClass, numClasses, outFileBase,
                         EMrep=NULL, distMetric=list(name="euclidean")){
  # deal with EMrep==NULL
  EMrep<-ifelse(is.null(EMrep),0,EMrep)
  # split off class number from row name
  classes <- as.numeric(sapply(strsplit(rownames(dataOrderedByClass),
                                        split="__class"),"[[",2))
  #dis<-stats::dist(dataOrderedByClass) # get distance matrix between reads
  dis<-getDistMatrix(dataOrderedByClass,distMetric)
  sil<-cluster::silhouette(classes,dis) # calculate silhouette
  # make data.frame with silhouette stats
  df<-data.frame(regionName=outFileBase, numClasses=numClasses, EMrep=EMrep,
                 silhouetteWidthMean=mean(sil[, 3],na.rm=T),
                 silhouetteWidthSD=stats::sd(sil[, 3],na.rm=T),
                 stringsAsFactors=F)

  classTable <- table(paste0("class",classes))

  # Add number of reads per class
  df[,paste0("class",1:numClasses,"_reads")]<-NA
  df[,paste0(names(classTable),"_reads")]<-classTable
  # Add average silhouette width per class
  df[,paste0("class",1:numClasses,"_silMean")]<-NA
  silWidthMean<-stats::aggregate(sil[,"sil_width"],by=list(sil[,"cluster"]),FUN=mean)
  colnames(silWidthMean)<-c("class","mean")
  df[,paste0("class",silWidthMean$class,"_silMean")]<-silWidthMean$mean
  #classMeans = stats::aggregate(dataMatrix, by = list(readClasses), FUN = mean)[-1]
  # Add silhouette width SD per class
  df[,paste0("class",1:numClasses,"_silSD")]<-NA
  silWidthSD<-stats::aggregate(sil[,3],by=list(sil[,1]),FUN=stats::sd)
  colnames(silWidthSD)<-c("class","sd")
  df[,paste0("class",silWidthSD$class,"_silSD")]<-silWidthSD$sd
  return(list(stats=df,plotObject=sil))
}


#' Plot class Means
#'
#' Create line plots for class means
#' @param classes A matrix of methylation or bincount values (classes x position) for each class
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param facet Plot mean profiles separately as a facet_wrap plot (default=TRUE)
#' @param title A title for the plot (default is "Class means")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object
#' @export
plotClassMeans<-function(classes,xRange=c(-250,250), facet=TRUE,
                         title="Class means",
                         myXlab="CpG/GpC position",featureLabel="TSS",
                         baseFontSize=12){
  # initialise variables
  position <- methFreq <- NULL
  numClasses<-nrow(classes)
  classMeans<-tidyr::gather(as.data.frame(classes),key="position",value="methFreq")
  classMeans$class<-as.factor(rep(paste0("class",1:numClasses),ncol(classes)))
  classMeans$position<-as.numeric(classMeans$position)


  p<-ggplot2::ggplot(classMeans,ggplot2::aes(x=position,y=1-methFreq,group=class)) +
    ggplot2::geom_line(ggplot2::aes(color=class))  +
    ggplot2::geom_point(ggplot2::aes(x=position,y=-0.07), size=0.5,
                        colour="grey80") +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("dSMF (1 - Methylation frequency)") +
    ggplot2::xlim(xRange) +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position="right",legend.box = "vertical",
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width=ggplot2::unit(0.3,"cm"))
  # add line for TSS
  p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=-0.07, ymax=1),
                               color="grey80") +
    ggplot2::annotate(geom="text", x=1,y=0.01,
                      label=featureLabel,color="grey20")
  if (facet==TRUE) {
    p<-p+ggplot2::facet_wrap(~class,nrow=nrow(classes))
  }
  return(p)
}




#' Plot all class means from multiple repeats
#'
#' Create line plots for class means
#' @param allClassMeans A table of positions, fraction methylation, classes for replicate EM runs
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param facet Plot mean profiles separately as a facet_wrap plot (default=TRUE)
#' @param title A title for the plot (default is "Class means of replicate EM runs")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object
#' @export
plotAllClassMeans<-function(allClassMeans,xRange=c(-250,250), facet=TRUE,
                            title="Class means of replicate EM runs",
                            myXlab="CpG/GpC position",featureLabel="TSS",
                            baseFontSize=12){
  # initialise variables
  position <- methFreq <- NULL
  numClasses<-max(as.integer(allClassMeans$class))
  allClassMeans$class<-as.factor(paste0("class",allClassMeans$class))
  allClassMeans$replicate<-as.factor(allClassMeans$replicate)
  allClassMeans$position<-as.numeric(allClassMeans$position)


  p<-ggplot2::ggplot(allClassMeans,ggplot2::aes(x=position, y=1-methFreq,
                                             group=replicate,
                                             colour=replicate)) +
    ggplot2::geom_line(ggplot2::aes(color=replicate),alpha=0.7)  +
    ggplot2::geom_point(ggplot2::aes(x=position,y=-0.07), size=0.5,
                        colour="grey80") +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("dSMF (1 - Methylation frequency)") +
    ggplot2::xlim(xRange) +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position="right",legend.box = "vertical",
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width=ggplot2::unit(0.3,"cm"))
  # add line for TSS
  p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=-0.07,ymax=1),
                               color="grey80") +
    ggplot2::annotate(geom="text", x=1,y=0.01,
                      label=featureLabel,color="grey20")
  if (facet==TRUE) {
    p<-p+ggplot2::facet_wrap(~class,nrow=numClasses)
  }
  return(p)
}



#' Plot smoothed class means from multiple repeats
#'
#' Plot loess-smoothed class means from multiple repeats of EM.
#' @param allClassMeans A long data frame of methylation or bincount values with columns for position, methylation frequency (methFreq), class and replicate.
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param facet Plot mean profiles separately as a facet_wrap plot (default=FALSE)
#' @param title A title for the plot (default is "Class means")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object
#' @export
plotSmoothedClassMeans<-function(allClassMeans, xRange=c(-250,250), facet=FALSE,
                                 title="Smoothed class means",
                                 myXlab="CpG/GpC position", featureLabel="TSS",
                                 baseFontSize=12){
  #initialise variables
  position <- methFreq <- NULL
  #numClasses<-max(as.integer(allClassMeans$class))
  allClassMeans$class<-as.factor(paste0("class",allClassMeans$class))
  allClassMeans$position<-as.numeric(allClassMeans$position)

  p<-ggplot2::ggplot(allClassMeans,ggplot2::aes(x=position,
                                                y=1-methFreq,
                                                group=class)) +
    ggplot2::geom_smooth(method="loess",span=0.5,
                         ggplot2::aes(color=class, fill=class)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("dSMF (1 - Methylation frequency)") +
    ggplot2::xlim(xRange) +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position="right",legend.box = "vertical")
  # add individual points if there are less than 100 positions (to distinguish
  # single molecule matrices from multigene matrices)
  if(length(unique(allClassMeans$position))<100){
    p<-p+ggplot2::geom_point(ggplot2::aes(colour=class),alpha=0.1)
  }
  # add line for TSS
  p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=0, ymax=1),
                               color="grey80") +
    ggplot2::annotate(geom="text", x=1, y=0.01,
                      label=featureLabel,color="grey20")
  return(p)
}


#' Plot smoothed class means from multiple repeats
#'
#' Plot loess-smoothed class means from multiple repeats of EM.
#' @param allClassMeans A long data frame of methylation or bincount values with columns for position, methylation frequency (methFreq), class and replicate.
#' @param regionGR A genomic range for the region that is being plotted
#' @param anchorPoint One of "middle" or "start": the position from which numbering starts for the position in allClassMatrix
#' @param txdb A TxDb object e.g. TxDb.Celegans.UCSC.ce11.refGene from which to get the gene regions
#' @param featureGR A genoimc range for a feature you wish to highlight (e.g. TSS). Will be overlaid on transcript track.
#' @param bigwigListFile A chacter string with the name of a file that contains a list of bigwigs to be plotted with the class means. First column of this file must be the file name of the bigwig with the full path. the second column is the dataset name (will be used on the plot so keep it short).
#' @param strandedBwFile A chacter string with the name of a file that contains a list of stranded bigwigs to be plotted with the class means. First column of this file must be the file name of the bigwig with the full path. The second column is the dataset name (will be used on the plot so keep it short), third column is "+" or "-" indicating the strand of the bigwig.
#' @return Returns a ggplot2 object
#' @export
plotClassMeansWithTracks<-function(allClassMeans, regionGR,
                                   anchorPoint="middle",
                                   txdb=NULL,
                                   featureGR=NULL,
                                   bigwigListFile=NULL,
                                   strandedBwFile=NULL){
  #initialise variables
  trackList<-grStrand<-grChr<-grGenome<-grStart<-grEnd<-NULL
  grStrand<-as.character(GenomicRanges::strand(regionGR))
  grChr<-as.character(GenomicRanges::seqnames(regionGR))
  grGenome<-GenomeInfoDb::genome(regionGR)[1]
  grStart<-GenomicRanges::start(regionGR)
  grEnd<-GenomicRanges::end(regionGR)

  numClasses=max(as.numeric(allClassMeans$class))

  itrack<-Gviz::IdeogramTrack(genome=grGenome, chromosome=grChr)

  genomeAxis <- Gviz::GenomeAxisTrack(name=grChr,labelPos="below")

  txTrack<-Gviz::GeneRegionTrack(txdb,name="Gene", chromosome=grChr,
                                 start=grStart, end=grEnd, genome=grGenome)

  gr<-allClassMeansToGR(allClassMeans, regionGR, anchorPoint="middle",
                              dSMFwin=1)
  GenomeInfoDb::genome(gr)<-grGenome
  dTrack<-Gviz::DataTrack(gr,name="dSMF Classes",
                          groups=rep(c(paste0("class",1:numClasses)),
                                     times=10),
                          type=c("a","p","confint"))

  trackList=c(itrack,genomeAxis)

  if(!is.null(featureGR)) {
    aTrack<-Gviz::AnnotationTrack(featureGR,
                                  shape="arrow",
                                  group=c(GenomicRanges::mcols(regionGR)$ID),
                                  groupAnnotation="group",
                                  just.group=ifelse(grStrand=="+",
                                                    "left","right"))
    oTrack<-Gviz::OverlayTrack(trackList=list(txTrack,aTrack))
    ht<-Gviz::HighlightTrack(trackList=list(oTrack,dTrack),
                           start=GenomicRanges::start(featureGR),
                           end=GenomicRanges::end(featureGR),
                           chromosome=grChr, genome=grGenome)

    trackList<-c(trackList,ht)
  } else {
    trackList=c(trackList,txTrack)
  }

  if(!is.null(strandedBwFile)){
    bigwigList<-utils::read.delim(strandedBwFile, header=FALSE,
                           stringsAsFactors=F)
    for (i in 1:nrow(bigwigList)){
      if(bigwigList[i,3]==grStrand){
        bw<-rtracklayer::import.bw(bigwigList[i,1], selection=regionGR)
        GenomeInfoDb::seqlevelsStyle(bw)<-"ucsc"
        bwTrack<-Gviz::DataTrack(bw, name=bigwigList[i,2], type="l")
        trackList<-c(trackList,bwTrack)
      }
    }
  }

  if(!is.null(bigwigListFile)){
    bigwigList<-utils::read.delim(bigwigListFile, header=FALSE,
                           stringsAsFactors=F)
    for (i in 1:nrow(bigwigList)){
      bw<-rtracklayer::import.bw(bigwigList[i,1], selection=regionGR)
      bwTrack<-Gviz::DataTrack(bw, name=bigwigList[i,2], type="hist")
      trackList<-c(trackList,bwTrack)
    }
  }

  p<-Gviz::plotTracks(trackList, from=grStart, to=grEnd, chromosome=grChr)
  return(p)
}


#' Convert allClassMeans table to GenomicRanges
#'
#' Convert allClassMeans table to GenomicRanges with absolute genomic
#' coordinates and with methylation score at each position for each class and
#' multiple replicates in the metadata columns (required for plotting with
#' Gviz)
#' @param allClassMeans A long data frame of methylation or bincount values with columns for position, methylation frequency (methFreq), class and replicate.
#' @param regionGR A genomic range for the region that is being plotted
#' @param anchorPoint One of "middle" or "start": the position from which numbering starts
#' @param dSMFwin Width of dSMF score window (default is 1bp)
#' @return Returns a GenomicRanges object with absolute genomic coordinates and
#' methylation score for different classes and replicates in mcols
#' @export
allClassMeansToGR<-function(allClassMeans, regionGR,
                            anchorPoint="middle",
                            dSMFwin=1){
  methFreq<-NULL
  numClasses=max(as.numeric(allClassMeans$class))
  allClassMeans$position<-as.numeric(allClassMeans$position)
  allClMeans<-tidyr::pivot_wider(allClassMeans,names_from=c(class,replicate),
                                 values_from=methFreq,names_prefix="class")
  grRelCoord<-GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(regionGR),
                             IRanges::IRanges(
                               start=as.numeric(allClMeans$position),
                               width=dSMFwin),
                             strand=GenomicRanges::strand(regionGR))
  # convert to absolute genomic coordinates
  if (anchorPoint=="middle") {
    starts<- GenomicRanges::start(regionGR) +
      GenomicRanges::width(regionGR)/2 + GenomicRanges::start(grRelCoord)
    ends<- GenomicRanges::start(regionGR) +
      GenomicRanges::width(regionGR)/2 + GenomicRanges::end(grRelCoord)
  } else if (anchorPoint=="start") {
    starts<- GenomicRanges::start(regionGR) - 1 +
      GenomicRanges::start(grRelCoord)
    ends<- GenomicRanges::start(regionGR) +
      GenomicRanges::end(grRelCoord)
  } else {
    print("anchorPoint must be one of 'middle' or 'start'")
  }
  gr<-GenomicRanges::GRanges(seqnames=GenomicRanges::seqnames(regionGR),
                                     IRanges::IRanges(start=starts,
                                       end=ends),
                                     strand=GenomicRanges::strand(regionGR))
  if(all(GenomicRanges::strand(gr)=="+")) {
    GenomicRanges::mcols(gr)<-1-allClMeans[,2:dim(allClMeans)[2]]
  } else if (all(GenomicRanges::strand(gr)=="-")) {
    GenomicRanges::mcols(gr)<-1-allClMeans[dim(allClMeans)[1]:1,2:dim(allClMeans)[2]]
  } else {
    "regionGR missing strand information"
  }
  return(gr)
}


#' Extract class info from dataOrderedByClass
#'
#' @param dataOrderedByClass A matrix of methylation frequency or bin counts for
#' indivudal reads at particular positions where the reads have been sorted by class and
#' the row names contain the read name and the class joined together: readName__classX
#' @param readNames A vector of read names by which to order the classes
#' @return Classification of reads ordered by readNames
#' @export
getReadClass<-function(dataOrderedByClass,readNames){
  readClasses<-as.numeric(sapply(strsplit(rownames(dataOrderedByClass),
                                           split="__class"),"[[",2))
  sortedReadNames<-sapply(strsplit(rownames(dataOrderedByClass),
                                              split="__class"),"[[",1)
  idx<-match(readNames,sortedReadNames)
  cl<-readClasses[idx]
  return(cl)
}


#' Get most frequently assigned class from repeat runs
#'
#' @param classVote A data.frame of reads with the class assignments from multiple EM repeats
#' @return A data.frame including the most frequently called class for each read
#' @export
getClassVote<-function(classVote){
  repCols<-names(classVote)[grep("rep",colnames(classVote))]
  if (length(repCols)>1) { # is there more than one repeat?
    classVote$topClass<-as.numeric(apply(classVote[,repCols],1,getMode))
    classVote$topClassFreq<-sapply(1:nrow(classVote),function(rowNum){
      topClassFreq<-sum(match(classVote[rowNum,repCols],classVote$topClass[rowNum]),
                      na.rm=T)/length(repCols)
      return(topClassFreq)})
  } else {
    classVote$topClass<--as.numeric(classVote[,repCols])
    classVote$topClassFreq<-1
  }
  return(classVote)
}



#' Plot single molecule and class means for individual repeats
#'
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced
#' @param outPath path to directory where plots will be saved
#' @param numClasses An integer indicating the number of classes to learn
#' @param EMrep An integer indicating which repeat of the EM this is
#' @param classMeans A matrix (classes x position) of the average methlyation profile of each class classes
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object
#' @export
plotEachRepeat<-function(dataOrderedByClass, outFileBase , outPath, numClasses,
                         EMrep, classMeans, xRange, myXlab, featureLabel,
                         baseFontSize){
  #makeDirs(path=outPath,dirNameList=paste0(c("classPlots",
  #                                         "classMeanPlots"),"/",outFileBase))
  # do single molecule plots of classes
  p<-plotClassesSingleMolecule(dataOrderedByClass, xRange,
                         title = outFileBase, myXlab=myXlab,
                         featureLabel=featureLabel, baseFontSize=baseFontSize)
  outPath<-gsub("\\/$","",outPath)
  ggplot2::ggsave(filename=paste0(outPath,
                                "/classReads_", outFileBase,"_K",
                                numClasses,"_r", EMrep, ".png"),
                plot=p, device="png", width=19, height=29, units="cm")

  #plot individual unsmoothed class means
  p<-plotClassMeans(classMeans,xRange=xRange, facet=TRUE,
                  title=paste(outFileBase, "Class means, repeat ", EMrep),
                  myXlab=myXlab, featureLabel=featureLabel,
                  baseFontSize=baseFontSize)

  ggplot2::ggsave(filename=paste0(outPath,
                                "/classMeans_", outFileBase,"_K",
                                numClasses, "_r",EMrep,".pdf"),
                plot=p, device="pdf", width=19, height=29, units="cm")
}



#' Plot histograms of consistency of classfiications over multiple repeats
#'
#' @param classVote A data.frame of reads with the class assignments from multiple EM repeats
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced
#' @param outPath path to directory where plots will be saved
#' @param numClasses An integer indicating the number of classes to learn
#' @return classVote data frame with "topClass" and "topClassFreq" columns added
#' @export
plotClassStability<-function(classVote,outFileBase,outPath,numClasses){
  #initialise varaibles
  topClassFreq <- topClass <- NULL
  p1<-ggplot2::ggplot(classVote,ggplot2::aes(x=topClassFreq, fill=topClass))+
    ggplot2::geom_histogram(binwidth=0.1)
  p2<-ggplot2::ggplot(classVote,ggplot2::aes(x=topClassFreq, fill=topClass))+
    ggplot2::geom_histogram(binwidth=0.1) + ggplot2::facet_wrap(.~topClass)

  p<-ggpubr::ggarrange(p1,p2,ncol=1,nrow=2)
  outPath<-gsub("\\/$","",outPath)
  ggplot2::ggsave(filename=paste0(outPath,"/classFreq_",
                                outFileBase,"_K",numClasses,".pdf"),
                plot=p,device="pdf",width=19,height=29,units="cm")
  return(classVote)
}



#' Draw silhouette plots and get silhouette statistics
#'
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param numClasses An integer indicating the number of classes to learn
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced
#' @param outPath path to directory where plots will be saved
#' @param EMrep An integer indicating which EM repeat this is
#' @param doIndividualPlots Produce individual plots for each repeat (default=F)
#' @param silStats data.frame with statistics about the silhouette plots (default=NULL)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return silStats data.frame with statistics about the silhouette plots
#' @export
getSilhouetteStats<-function(dataOrderedByClass, numClasses, outFileBase, outPath,
                             EMrep=NULL, doIndividualPlots=FALSE, silStats=NULL,
                             distMetric=list(name="euclidean")){
  silList<-silhouettePlot(dataOrderedByClass, numClasses, outFileBase,
                          EMrep, distMetric)

  if (!is.null(silStats)) {
    silStats<-rbind(silStats,silList$stats)
  } else {
    silStats<-silList$stats
  }

  # save silhouette for individual repeats if first round or if requested
  if(is.null(EMrep) | doIndividualPlots==TRUE) {
    #makeDirs(path=outPath,dirNameList=paste0("silPlts","/",outFileBase))
    repTxt<-ifelse(is.null(EMrep),"",paste0("_r",EMrep))
    outPath<-gsub("\\/$","",outPath)
    grDevices::pdf(paste0(outPath,"/sil_",
                          outFileBase,"_K",
                          numClasses, repTxt, ".pdf"),
                   paper="a4", height=11, width=8)
    graphics::plot(silList$plotObject,main=paste0("Silhouette: ",
    outFileBase," with ",numClasses," classes",repTxt))
    graphics::abline(v=silList$stats["silhouetteWidthMean"], col="black",lty=2)
    grDevices::dev.off()
  }
  return(silStats)
}


#' Plot reads and class means with final classification
#'
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param numClasses An integer indicating the number of classes to learn
#' @param allClassMeans A long data frame of methylation or bincount values with columns for position, methylation frequency (methFreq), class and replicate.
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced
#' @param outPath path to directory where plots will be saved
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return None
#' @export
plotFinalClasses<-function(dataOrderedByClass, numClasses, allClassMeans,
                           outFileBase,
                           outPath, xRange, myXlab, featureLabel, baseFontSize){
  # plot single molecules with final classes
  p<-plotClassesSingleMolecule(dataOrderedByClass, xRange=xRange,
                         title = outFileBase, myXlab=myXlab,
                         featureLabel=featureLabel, baseFontSize=12)
  outPath<-gsub("\\/$","",outPath)
  ggplot2::ggsave(filename=paste0(outPath,
                                "/finalClassReads_", outFileBase,"_K",
                                numClasses, ".png"),
                plot=p, device="png", width=19, height=29, units="cm")


  # plot all class means (faceted)
  p<-plotAllClassMeans(allClassMeans,xRange=xRange, facet=TRUE,
                              title=paste0(outFileBase, " class means of ",
                                           max(allClassMeans$replicate),
                                           " EM runs"),
                              myXlab="CpG/GpC position", featureLabel="TSS",
                              baseFontSize=12)

  ggplot2::ggsave(filename=paste0(outPath,"/allClassMeans_",
                                  outFileBase,"_K", numClasses,".pdf"),
                  plot=p, device="pdf", width=19, height=29, units="cm")


  # plots smoothed average class means +- StdErr
  EMrepeats<-max(allClassMeans$replicate)
  p<-plotSmoothedClassMeans(allClassMeans, xRange=xRange, facet=TRUE,
                          title=paste(outFileBase, " Class means, ",
                                      EMrepeats ," repeats"),
                          myXlab="CpG/GpC position", featureLabel="TSS",
                          baseFontSize=12)

  ggplot2::ggsave(filename=paste0(outPath,"/smClassMeans_",
                                outFileBase,"_K", numClasses,".pdf"),
                plot=p, device="pdf", width=25, height=12, units="cm")
}





#' Run several repeats of iterative EM
#'
#' Perform EM clustering on the same matrix several times. Necessary to check ¨
#' that classes are stable despite the random assignment of methylation
#' fractions and the random assignment of starting conditions.
#'
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param EMrepeats An integer indicating the number of times to repeat the clustering (default=10)
#' @param outPath A string with the path to the directory where the output should go
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced (default is "")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @param doIndividualPlots Produce individual plots for each repeat (default=F)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return  allClassMeans data.frame with columns: position, methFreq, class, replicate
#' @export
runEMrepeats<-function(dataMatrix, numClasses=3, convergenceError=1e-6,
                       maxIterations=100, EMrepeats=10, outPath=".",
                       xRange=c(-250,250), outFileBase="",
                       myXlab="CpG/GpC position", featureLabel="TSS",
                       baseFontSize=12, doIndividualPlots=FALSE,
                       distMetric=list(name="euclidean")){
  #initialise variables
  methFreq <- position <- NULL
  #remove trailing / from outPath
  outPath<-gsub("\\/$","",outPath)
  # make output directories
  #makeDirs(path=outPath,dirNameList=c("silPlts","dataOrdByClass","classPlots",
  #                                    "classMeanPlots"))

  previousClassMeans<-NULL # in the first round, use hclust to sort clusters
  classVote<-data.frame(read=row.names(dataMatrix),stringsAsFactors=F)
  silStats<-NULL

  # check if repeats necessary for this matrix
  #fractionIdx<-dataMatrix > 0 & dataMatrix < 1
  dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
  stopifnot(isMatrixValid(dataMatrix, NAsValid=FALSE))

  for (EMrep in 1:EMrepeats){
    # do classifiction
    emClass<-runEM(dataMatrix=dataMatrix, numClasses=numClasses,
                   convergenceError=convergenceError,
                   maxIterations=maxIterations)

    # order dataMatrix by class
    orderedData<-classifyAndSortReads(dataMatrix=dataMatrix,
                                      posteriorProb=emClass$posteriorProb,
                                      previousClassMeans=previousClassMeans)
    dataOrderedByClass<-orderedData$data
    classMeans<-orderedData$classMeans

    classVote[,paste0("rep",EMrep)]<-getReadClass(dataOrderedByClass,
                                                  classVote$read)

    #previousClassMeans<-setPreviousClassMeans(classMeans, previousClassMeans, allClassMeans)
    # store classMeans from first round as previousClassMeans to have consitent order
    # store classMeans from all rounds in allClassMeans in order to return at end of function
    if(is.null(previousClassMeans)) {
      previousClassMeans<-classMeans
      allClassMeans<-tidyr::gather(classMeans,key="position",value="methFreq")
      allClassMeans$class<-rep(rownames(classMeans),ncol(classMeans))
      allClassMeans$replicate<-EMrep
    } else {
      tmpMeans<-tidyr::gather(classMeans,key="position",value="methFreq")
      tmpMeans$class<-rep(rownames(classMeans),ncol(classMeans))
      tmpMeans$replicate<-EMrep
      allClassMeans<-rbind(allClassMeans,tmpMeans)
    }

    # plot the results for individual repeats if first round and/or if requested
    if(doIndividualPlots==TRUE) { #if want 1 repeat can add:  | EMrep==1
      print("plotting individual EMruns")
      plotEachRepeat(dataOrderedByClass, outFileBase, outPath, numClasses,
                     EMrep,
                     classMeans, xRange, myXlab, featureLabel, baseFontSize)

    }
    print("getting silhouette statistics")
    # do silhouette plot and save silhouette stats
    silStats<-getSilhouetteStats(dataOrderedByClass, numClasses, outFileBase,
                                 outPath, EMrep, doIndividualPlots, silStats,
                                 distMetric)
  }

  #print("plotting class stability")
  # plots barplots of consistency of classfiications over multiple repeats
  classVote<-getClassVote(classVote)
  classVote$topClass<-factor(classVote$topClass)
  plotClassStability(classVote,outFileBase,outPath,numClasses)

  # save data with most frequent class call.
  #print("saving data with most frequent class call")
  idx<-match(row.names(dataMatrix),classVote$read)
  row.names(dataMatrix)<-paste0(rownames(dataMatrix),"__class",
                                classVote$topClass[idx])
  dataOrderedByClassRep<-dataMatrix[order(classVote$topClass[idx]),]

  saveRDS(dataOrderedByClassRep, file=paste0(outPath, "/", outFileBase, "_K",
                                             numClasses, ".rds"))
  print("plotting final classes")
  plotFinalClasses(dataOrderedByClassRep, numClasses, allClassMeans,
                   outFileBase,
                   outPath, xRange, myXlab, featureLabel, baseFontSize)

  #print("plotting silhouette statistics")
  # do silhouette plot and save silhouette stats
  silStats<-getSilhouetteStats(dataOrderedByClassRep, numClasses, outFileBase,
                               outPath, EMrep=NULL, silStats=silStats,
                               distMetric=distMetric)

  #calculate elbow and gap statistic
  #print("calculating gap statistic")
  readClasses <- paste0("class",sapply(strsplit(rownames(dataOrderedByClassRep),
                                 split="__class"),"[[",2))
  silStats$elbowWSS<-withinClusterSS(dataOrderedByClassRep, readClasses,
                                     distMetric)

  #print("saving silouhette stats")
  utils::write.csv(silStats, file=paste0(outPath,"/silStats_",
                                  outFileBase,"_K",numClasses,".csv"),
                   row.names=F)
  return(allClassMeans)
}



#' Run EM to predict a range of class numbers
#'
#' Perform EM clustering on the same matrix with a range of different class numbers.
#'
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param k_range A vector indicating different numbers of classes to learn
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param EMrepeats An integer indicating the number of times to repeat the clustering (default=10)
#' @param outPath A string with the path to the directory where the output should go
#' @param xRange A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced (default is "")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @param doIndividualPlots Produce individual plots for each repeat (default=F)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return allClassMeans list of different class numbers each containing a data.frame with columns: position, methFreq, class, replicate
#' @export
runEMrangeClassNum<-function(dataMatrix, k_range=2:8, convergenceError=1e-6,
                             maxIterations=100, EMrepeats=10, outPath=".",
                             xRange=c(-250,250), outFileBase="",
                             myXlab="CpG/GpC position", featureLabel="TSS",
                             baseFontSize=12, doIndividualPlots=TRUE,
                             distMetric=list(name="euclidean")) {
  #stopifnot(isMatrixValid(dataMatrix))
  allClassMeans<-list()
  for (numClasses in k_range) {
    print(paste("numClasses:",numClasses))
    allClassMeans[[numClasses]]<-runEMrepeats(dataMatrix, numClasses,
                                            convergenceError, maxIterations,
                                            EMrepeats, outPath, xRange,
                                            outFileBase, myXlab, featureLabel,
                                            baseFontSize, doIndividualPlots,
                                            distMetric)
  }
  return(allClassMeans)
}



#' Run several repeats of iterative EM for randomised data
#'
#' Perform EM clustering on the same matrix several times. Necessary to check that classes are stable despite the random assignment of methylation fractions.
#'
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param EMrepeats An integer indicating the number of times to repeat the clustering (default=10)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return  Numeric. Total within cluster sum of squares.
#' @export
runEMrepeats_withinSS<-function(dataMatrix, numClasses=3, convergenceError=1e-6,
                                maxIterations=100, EMrepeats=10,
                                distMetric=list(name="euclidean")){
  previousClassMeans<-NULL # in the first round, use hclust to sort clusters
  classVote<-data.frame(read=row.names(dataMatrix),stringsAsFactors=F)
  for (EMrep in 1:EMrepeats) {
    # do classifiction
    emClass<-runEM(dataMatrix=dataMatrix, numClasses=numClasses,
                   convergenceError=convergenceError,
                   maxIterations=maxIterations, printProgress=FALSE)

    # order dataMatrix by class
    orderedData<-classifyAndSortReads(dataMatrix=dataMatrix,
                                      posteriorProb=emClass$posteriorProb,
                                      previousClassMeans=previousClassMeans)
    dataOrderedByClass<-orderedData$data
    classMeans<-orderedData$classMeans
    classVote[,paste0("rep",EMrep)]<-getReadClass(dataOrderedByClass,
                                                  classVote$read)

    # store classMeans from first round as previousClassMeans to have consitent order
    if(is.null(previousClassMeans)) {
      previousClassMeans<-classMeans
    }
  }

  # plots histograms of consistency of classfiications over multiple repeats
  classVote<-getClassVote(classVote)
  classVote$topClass<-factor(classVote$topClass)

  # calculate total within cluster sum of squares
  idx<-match(row.names(dataMatrix),classVote$read)
  classes<-classVote$topClass[idx]
  withinSS<-withinClusterSS(dataMatrix,classes,distMetric)
  return(withinSS)
}



#' Plot clustering metrics
#'
#' Diagnostic plots for choosing the optimal number of clusters. Based on the following:
#' https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/. Using three metrics:
#' 1) average silhouette width. Bigger is better
#' 2) Within cluster sum of squares (elbowWSS) - smaller is better, look for where slope no longer changes much ("elbow" in graph)
#' 3) Gap statistic which compares how much smaller the sample WSS is from a random distribution of WSS generated from randomised matrices. Bigger is better. Look for kink in graph.
#' @param dataMatrix A matrix of methylation or bincount values (reads x position)
#' @param k_range A vector indicating different numbers of classes to learn
#' @param maxB  The maximum number of randomisations to perform
#' @param convergenceError An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param outPath A string with the path to the directory where the output should go
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced (default is "")
#' @param EMrep An integer indicating which EM repeat this is
#' @param nThreads Number of threads to use for generating background distribution (default is 1)
#' @param setSeed Logical value to determine if seed should be set for randomisation (default is FALSE)
#' @param distMetric A list with the name of the distance metric and any
#' parameters it might require
#' @return None
#' @export
plotClusteringMetrics<-function(dataMatrix, k_range=2:8, maxB=100,
                                convergenceError=1e-6, maxIterations=100,
                                outPath=".", outFileBase="", EMrep=NULL,
                                nThreads=1, setSeed=FALSE,
                                distMetric=list(name="euclidean")){
  # initialise variables
  meanSilhouetteWidth <- elbowWSS <- gap <- position <- NULL
  classMean <- classNumber <- NULL
  outPath<-gsub("\\/$","",outPath)
  EMrep<-ifelse(is.null(EMrep),0,EMrep)

  clusterStats<-data.frame(numClasses=k_range, meanSilhouetteWidth=NA,
                           elbowWSS=NA, gap=NA,
                           stringsAsFactors=F)
  # add columns for individual class silhouette widths
  individClassSil<-data.frame(matrix(data=NA,nrow=length(k_range),
                                     ncol=max(k_range)),stringsAsFactors=F)
  colnames(individClassSil)<-paste0("class",1:max(k_range),"_silMean")
  clusterStats<-cbind(clusterStats,individClassSil)

  print("generating randomised matrices as null distribution")
  randomisedMatrixStatsFile<-paste0(outPath,"/randMatStats_",
                                    outFileBase,".csv")
  if (file.exists(randomisedMatrixStatsFile)) {
    randomWSS<-utils::read.csv(randomisedMatrixStatsFile, stringsAsFactors=F)
  } else {
    randomWSS<-clusterRandomMatrices(dataMatrix, k_range, maxB,
                                     convergenceError, maxIterations,
                                     nThreads, setSeed, distMetric)
    utils::write.csv(randomWSS, randomisedMatrixStatsFile, row.names=F)
  }

  print("calculating summary statistics for clustering")
  for (numClasses in k_range) {
    print(numClasses)
    silStats<-utils::read.csv(paste0(outPath,"/silStats_",
                            outFileBase,"_K",numClasses,".csv"),stringsAsFactors=F)

    nc<-which(clusterStats$numClasses==numClasses)
    # add individual class mean silhouette widths
    classSilMeanCols<-silStats[silStats$EMrep==EMrep,
                               grep("_silMean",colnames(silStats))]
    classSilMeanCols[is.na(classSilMeanCols)]<-0
    clusterStats[nc,colnames(classSilMeanCols)]<-classSilMeanCols

    # add overall silhouette mean width
    clusterStats$meanSilhouetteWidth[nc]<-
      silStats$silhouetteWidthMean[silStats$EMrep==EMrep]

    clusterStats$elbowWSS[nc]<-mean(silStats$elbowWSS)
    clusterStats$gap[nc]<-
      randomWSS[randomWSS$numClasses==numClasses,"sumSq"]-clusterStats$elbowWSS[nc]
  }

  classMeanCols<-colnames(clusterStats)[grep("_silMean",colnames(clusterStats))]
  long<-tidyr::gather(clusterStats,key="classNumber",
                      value="classMean", classMeanCols)
  long$classNumber<-gsub("_silMean","",long$classNumber)

  print("plotting clustering statistics")
  p1<-ggplot2::ggplot(long,ggplot2::aes(x=numClasses,
                                                y=meanSilhouetteWidth)) +
    ggplot2::geom_line(ggplot2::aes(x=numClasses,y=meanSilhouetteWidth)) +
    ggplot2::geom_point(ggplot2::aes(x=numClasses,y=classMean,
                                     colour=classNumber), alpha=0.5) +
    ggplot2::geom_hline(yintercept=0, colour="red", size=0.2)+
    ggplot2::ggtitle(paste("Silhouette width", outFileBase))
  #p2<-ggplot2::ggplot(clusterStats,ggplot2::aes(x=numClasses,y=elbowWSS)) +
  #  ggplot2::geom_line() + ggplot2::geom_point() +
  #  ggplot2::ggtitle(paste("Within class Sum of Squares", outFileBase))
  p3<-ggplot2::ggplot(clusterStats,ggplot2::aes(x=numClasses,y=gap)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ggtitle(paste("gap statistic", outFileBase))
  #p<-ggpubr::ggarrange(p1,p2,p3,nrow=3,ncol=1)
  p<-ggpubr::ggarrange(p1,p3,nrow=2,ncol=1)
  ggplot2::ggsave(paste0(outPath,"/clustStats_",
                         outFileBase,".pdf"),
                  plot=p, device="pdf", width=19,height=29,units="cm")
  return("Clustering metrics plotted successfully")
}




