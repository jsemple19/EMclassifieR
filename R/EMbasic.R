#' Expectation Maximisation algorithm
#'
#' @param classes A matrix containing the classes to be optimised. c[i,j] is the expected methylation value or bincount value of class i at position j.
#' @param priorProb A vector defining the prior probabilities of each class.
#' @param data  A matrix containing the sample. data [i,j] is the methylation value or bincount value of sample i at position j.
#' @return A list of three items:
#'             1) classes: a matrix with classes being optimised (class x position)
#'             2) priorProb: a vector with the prior probabilities of each class
#'             3) posteriorProb: a matrix of probabilites of each sample belonging
#'             to a particular class (samples x class)
#' @export
em_basic <- function(classes,priorProb,data) {
  # To deal with fractions in methlyation data (when opposite strand do not agree on
  # methylation status), the methylation status will be applied at random with
  # probability given by the fraction.
  fractionIdx<-data > 0 & data < 1
  binaryData<-data
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx),size=1,prob=data[fractionIdx])
  }

  numClasses=dim(classes)[1]         # Number of classes
  numSamples=dim(binaryData)[1]            # Number of samples
  l=matrix(nrow=numSamples, ncol=numClasses)  # log linumClasseselihood
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability by which each class occurs

  # E step
  for(i in 1:numSamples) {
    for (j in 1:numClasses) {
      # Bernoulli (size = 1): simply multiplying by class probability
      l[i,j] = sum(stats::dbinom(x = binaryData[i,], size = 1, prob = classes[j,], log = T))
      # not sure about implementing the beta distribution
      #l[i,j] = sum(dbeta(x = data[i,], shape1 = classes[j,], shape2 = classes[j,], log = T))
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



############################### EM HELPERS #########################################


#' Check if p value has converged
#'
#' @param p_diff The difference in likelihood betwween this round of optimisation and the previous round.
#' @param p_diff_prev The difference in likelihood between the previous round of optimisation and the one before that.
#' @param error The maximum tolerated error
#' @return A TRUE or FALSE value indicating if the difference between successive rounds of optimisation is below the error threshold (i.e. converged)
p_converged <- function(p_diff, p_diff_prev, error) {
  return(abs(p_diff-p_diff_prev) < error)
}


#' Get absolute value of difference between two p values
#'
#' @param p Current probability
#' @param p_prev Probability from last round of optimisation
#' @return The absolute value of the change in probability
p_difference <- function(p, p_prev) {
  return(sum(abs(p-p_prev)))
}



#' Run EM iteratively to convergence
#'
#' @param data A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergence_error An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @return  list of three items:
#'             1) classes: a matrix with optimised classes (class x position)
#'             2) priorProb: a vector with the prior probabilities of each class
#'             3) posteriorProb: a matrix of probabilites of each sample belonging
#'             to a particular class (samples x class)
#' @export
runEM<-function(data, numClasses, convergence_error=1e-6, maxIterations=100) {
  ################################# RUN EM #########################################
  #
  # Complete algorithm for partitioning with random seeds

  # To deal with fractions in methlyation data (when opposite strand do not agree on
  # methylation status), the methylation status will be applied at random with
  # probability given by the fraction.
  fractionIdx<-data > 0 & data < 1
  binaryData<-data
  if (sum(fractionIdx)>0){
    binaryData[fractionIdx]<-stats::rbinom(n=sum(fractionIdx),size=1,prob=data[fractionIdx])
  }

  numSamples=dim(binaryData)[1]            # Number of samples
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability each sample (read) belongs to a particular class

  # set.seed(12)    # same random probabilities each time
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


  e = convergence_error * length(posteriorProb) # scale convergence error to size of p

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
  print(p_converged(p_diff, p_diff_prev, e) & i < maxIterations)
  while (!p_converged(p_diff, p_diff_prev, e) & i < maxIterations) {
    print(paste("i:", i, " , sum:", p_diff, " , step diff:", abs(p_diff-p_diff_prev)))
    classes[classes>1]=1 # dbern error when probability > 1
    p_prev = posteriorProb
    p_diff_prev = p_diff

    results<-em_basic(classes,priorProb,binaryData)
    classes=results$classes
    posteriorProb=results$posteriorProb
    priorProb=results$priorProb
    p_diff = p_difference(posteriorProb, p_prev)

    i = i+1
    #print(posteriorProb)
    #print(p_diff)
    #print(p_diff_prev)
  }
  print(paste("converged:", p_converged(p_diff, p_diff_prev, e) & i < maxIterations, "iterations: ",i))
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
#' If the means of previous classes was not provided, hclust is used to cluste classes by their similarity.
#' @param data A matrix of methylation or bincount values (reads x position)
#' @param posteriorProb posteriorProb: a matrix of probabilites of each sample belonging to a particular class (samples x class)
#' @param previousClassMeans A matrix of the class means from a previous round of clustering
#' @return Returns a matrix with the reads classified (__classX is appended to the read name), and the classes are sorted.
#' @export
classifyAndSortReads<-function(data,posteriorProb,previousClassMeans=NULL) {
  ###################### POST-EM DATA SORTING ###########################################
  #
  # assign classes to reads according to the highest class probability
  numClasses=ncol(posteriorProb)
  readClasses = apply(posteriorProb, 1, which.max)
  #readsTable = table(readClasses)
  classMeans = stats::aggregate(data, by = list(readClasses), FUN = mean)[-1]

  if (!is.null(previousClassMeans)) {
    print("orderByPreviousClusters")
    # order the classes by comparing the class means to the previous clustering
    classOrder = order_by_prev_cluster(numClasses, classMeans, previousClassMeans)
    classMeans = classMeans[classOrder, ]
    rownames(classMeans)<-c(1:numClasses)[c(1:numClasses) %in%
                                            unique(readClasses)]
  } else {
    print("orderByHClust")
    # order the classes by similarity (class means)
    classOrder = stats::hclust(stats::dist(classMeans))$order
    classMeans = classMeans[classOrder, ]
    rownames(classMeans)<-c(1:numClasses)[c(1:numClasses) %in%
                                            unique(readClasses)]
  }
  print(classOrder)
  # change name of classes to match the order
  readClasses<-factor(readClasses)
  # remap the values in the vector
  rc<-plyr::mapvalues(readClasses,
                      from = levels(readClasses)[classOrder],
                      to = levels(readClasses))
  # change the level order
  rc1<-factor(rc,levels=c(1:numClasses))
  rownames(data)<-paste(rownames(data),rc1,sep="__class")
  ord = order(match(rc1,1:numClasses))
  dataOrderedByClass = data[ord,]

  return(list(data=dataOrderedByClass, classMeans=classMeans))
}



#' Plot reads by class for a single gene
#'
#' Create a single molecule plot of the reads sorted by class.
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param xlim A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param title A title for the plot (default is "Reads by classes")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object of a single molecule plot sorted by classes
#' @export
plotClassesSingleGene<-function(dataOrderedByClass,
                                xlim=c(-250,250), title="Reads by classes",
                                myXlab="CpG/GpC position",
                                featureLabel="TSS", baseFontSize=12){
  readClasses <- sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",2)
  classOrder <- unique(readClasses)
  readNames<-sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",1)
  readsTable <- table(readClasses)
  #classMeans = stats::aggregate(data, by = list(readClasses), FUN = mean)[-1]
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
  d$position<-as.numeric(d$position)
  p<-ggplot2::ggplot(d,ggplot2::aes(x=position,y=molecules)) +
    ggplot2::geom_tile(ggplot2::aes(width=4,fill=methylation),alpha=0.8) +
    ggplot2::scale_fill_gradient(low="blue", high="red", na.value="transparent",
                                 breaks=c(0,1), labels=c("protected","accessible"),
                                 limits=c(0,1), name="dSMF\n\n") +
    #ggplot2::scale_fill_manual(values=c("0"="black","1"="grey80"),na.translate=F,na.value="white", labels=c("protected","accessible"),name="dSMF") +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.position="right",legend.box = "vertical",
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width=ggplot2::unit(0.3,"cm")) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("Single molecules") +
    ggplot2::xlim(xlim[1],xlim[2]+10)
    # add line for TSS
    p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=0,
                                              ymax=length(reads) +
                                                max(3, 0.04*length(reads))),
                                              color="grey80")+
      ggplot2::annotate(geom="text", x=1,
                        y=-max(2,0.03*length(reads)),
                        label=featureLabel,color="grey20")
    # add lines separating classes
    p<-p+ggplot2::geom_hline(yintercept=classBorders,colour="grey80")
    # add color bar for classes
    p<-p+ggplot2::geom_segment(data=df1,mapping=ggplot2::aes(x=xlim[2]+10,y=readNumber-0.5,xend=xlim[2]+10,yend=readNumber+0.5,colour=Class),size=5)+
      ggplot2::geom_vline(xintercept=xlim[2],colour="grey80")
    return(p)
}



#' Plot silhouette plot to evaluate classification
#'
#' Create a silhouette plot to evaluate the quality of classification, and return some basic parameters about the classification
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param numClasses An integer with the number of classes learned
#' @param silhouetteDir String denoting the directory in which to save silhouette plots
#' @param regionName String with the name of the region for which reads are being classified.
#' @return Creates a silhoutte plot in a pdf file and returns a dataframe with mean and SD of silhouette width overall, and per class, as well as number of reads per class
#' @export
silhouettePlot<-function(dataOrderedByClass, numClasses, silhouetteDir,
                                 regionName){
  # split off class number from row name
  classes <- as.numeric(sapply(strsplit(rownames(dataOrderedByClass),
                                        split="__class"),"[[",2))
  dis<-stats::dist(dataOrderedByClass) # get distance matrix between reads
  sil<-cluster::silhouette(classes,dis) # caluculate silhouette
  # makd data.frame with silhouette stats
  df<-data.frame(regionName=regionName,numClasses=numClasses,
                 silhouetteWidthMean=mean(sil[, 3],na.rm=T),
                 silhouetteWidthSD=stats::sd(sil[, 3],na.rm=T),stringsAsFactors=F)
  grDevices::pdf(paste(silhouetteDir,"/", regionName, "_silhouette_K", numClasses,
            ".pdf", sep=""),paper="a4",height=11,width=8)
  graphics::plot(sil)
  graphics::abline(v=df$silhouetteWidthMean, col="black",lty=2)
  grDevices::dev.off()

  classTable <- table(paste0("class",classes))

  # Add number of reads per class
  df[,paste0("class",1:numClasses,"_reads")]<-NA
  df[,paste0(names(classTable),"_reads")]<-classTable
  # Add average silhouette width per class
  df[,paste0("class",1:numClasses,"_silMean")]<-NA
  silWidthMean<-stats::aggregate(sil[,3],by=list(sil[,1]),FUN=mean)
  colnames(silWidthMean)<-c("class","mean")
  df[,paste0("class",silWidthMean$class,"_silMean")]<-silWidthMean$mean
  #classMeans = stats::aggregate(data, by = list(readClasses), FUN = mean)[-1]
  # Add silhouette width SD per class
  df[,paste0("class",1:numClasses,"_silSD")]<-NA
  silWidthSD<-stats::aggregate(sil[,3],by=list(sil[,1]),FUN=stats::sd)
  colnames(silWidthSD)<-c("class","sd")
  df[,paste0("class",silWidthSD$class,"_silSD")]<-silWidthSD$sd
  return(df)
}


#' Plot class Means
#'
#' Create a single molecule plot of the reads sorted by class.
#' @param classes A matrix of methylation or bincount values (classes x position) for each class
#' @param xlim A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param title A title for the plot (default is "Class means")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param overplot Plot mean profiles separately as a facet_wrap plot (default=TRUE).
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return Returns a ggplot2 object of class means
#' @export
plotClassMeans<-function(classes,xlim=c(-250,250), facet=TRUE, title="Class means",
                         myXlab="CpG/GpC position",featureLabel="TSS",
                         baseFontSize=12){
  numClasses<-nrow(classes)
  classMeans<-tidyr::gather(as.data.frame(classes),key="position",value="methFreq")
  classMeans$class<-as.factor(rep(paste0("class",1:numClasses),ncol(classes)))
  classMeans$position<-as.numeric(classMeans$position)


  p<-ggplot2::ggplot(classMeans,ggplot2::aes(x=position,y=1-methFreq,group=class)) +
    ggplot2::geom_line(ggplot2::aes(color=class))  +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(myXlab) +
    ggplot2::ylab("dSMF (1 - Methylation frequency)") +
    ggplot2::xlim(xlim[1],xlim[2]+10) +
    ggplot2::theme_light(base_size=baseFontSize) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.position="right",legend.box = "vertical",
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width=ggplot2::unit(0.3,"cm"))
  # add line for TSS
  p<-p+ggplot2::geom_linerange(ggplot2::aes(x=1, y=NULL, ymin=0,ymax=1),
                               color="grey80") +
    ggplot2::annotate(geom="text", x=1,y=0.01,
                      label=featureLabel,color="grey20")
  if (facet==TRUE) {
    p<-p+ggplot2::facet_wrap(~class,nrow=nrow(classes))
  }
  return(p)
}






#' Run several repeats of iterative EM
#'
#' Perform EM clustering on the same matrix several times. Necessary to check that classes are stable despite the random assignment of methylation fractions.
#'
#' @param data A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergence_error An float indicating the convergence threshold for stopping iteration
#' @param maxIterations An integer indicating the max number of iterations to perform even if the algorithm has not converged
#' @param repeats An integer indicating the number of times to repeat the clustering (default=10)
#' @param outPath A string with the path to the directory where the output should go
#' @param xlim A vector of the first and last coordinates of the region to plot (default is c(-250,250))
#' @param outFileBase A string that will be used in the filenames and titles of the plots produced (default is "")
#' @param myXlab  A label for the x axis (default is "CpG/GpC position")
#' @param featureLabel A label for a feature you want to plot, such as the position of the TSS (default="TSS")
#' @param baseFontSize The base font for the plotting theme (default=12 works well for 4x plots per A4 page)
#' @return  None (value 0)
#' @export
runEMrepeats<-function(data, numClasses=3, convergence_error=1e-6, maxIterations=100,
                       repeats=10, outPath=".", xlim=c(-250,250), outFileBase = "",
                       myXlab="CpG/GpC position", featureLabel="TSS",
                       baseFontSize=12){
  # make output directories
  makeDirs(path=outPath,dirNameList=c("silhouettePlots","dataOrderedByClass","classPlots",
                                      "classMeanPlots"))
  previousClassMeans=NULL
  for (rep in 1:repeats) {
    # do classifiction
    emClass<-runEM(data=data, numClasses=numClasses, convergence_error=convergence_error,
                   maxIterations=maxIterations)

    # order data by class
    orderedData<-classifyAndSortReads(data=data, posteriorProb=emClass$posteriorProb,
                                      previousClassMeans=previousClassMeans)
    dataOrderedByClass<-orderedData$data
    classMeans<-orderedData$classMeans

    if(is.null(previousClassMeans)) {
      previousClassMeans<-classMeans
    }

    saveRDS(dataOrderedByClass, file=paste0(outPath, "/dataOrderedByClass/",
                                            outFileBase,"_rep",rep,".pdf"))

    # do single molecule plots of classes
    p<-plotClassesSingleGene(dataOrderedByClass=dataOrderedByClass, xlim=xlim,
                             title = outFileBase, myXlab=myXlab,
                             featureLabel=featureLabel, baseFontSize=12)


    ggplot2::ggsave(filename=paste0(outPath,"/classPlots/classifiedReads_",
                                    outFileBase,"_rep", rep, ".pdf"),
                    plot=p, device="pdf", width=19, height=29, units="cm")

    # do line plots of class mean values
    p<-plotClassMeans(classMeans,xlim=xlim, facet=FALSE,
                      title=paste(outFileBase, "Class means, repeat ", rep),
                      myXlab=myXlab, featureLabel=featureLabel,
                      baseFontSize=12)

    ggplot2::ggsave(filename=paste0(outPath,"/classMeanPlots/classMeans_",
                                    outFileBase, "_rep", rep, ".pdf"),
                    plot=p, device="pdf", width=29, height=19, units="cm")

    p<-plotClassMeans(classMeans,xlim=xlim, facet=TRUE,
                      title=paste(outFileBase, "Class means, repeat ", rep),
                      myXlab=myXlab, featureLabel=featureLabel,
                      baseFontSize=12)

    ggplot2::ggsave(filename=paste0(outPath,"/classMeanPlots/classMeans_facet_",
                                    outFileBase,"_rep",rep,".pdf"),
                    plot=p, device="pdf", width=19, height=29, units="cm")



    # do silhouette plot and silhouette data
    df<-silhouettePlot(dataOrderedByClass, numClasses,
                       silhouetteDir=paste0(outPath,"/silhouettePlots"),
                       regionName=paste0(outFileBase))
    write.csv(df, file=paste0(outPath,"/silhouettePlots/silhouetteStats_",
                              outFileBase,"_rep",rep,".csv"))
  }
  return(0)
}


#runDifferentClassSizes()

#runManyGenesTogether()
