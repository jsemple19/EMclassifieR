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
#' @param p_diff
#' @param p_diff_prev
#' @param error
#' @return A TRUE or FALSE value
p_converged <- function(p_diff, p_diff_prev, error) {
  return(abs(p_diff-p_diff_prev) < error)
}

#' Get absolute value of difference between two p values
#'
#' @param p
#' @param p_prev
#' @param error
#' @return The absolute value of the difference
p_difference <- function(p, p_prev) {
  return(sum(abs(p-p_prev)))
}



#' Run EM iteratively
#'
#' @param data A matrix of methylation or bincount values (reads x position)
#' @param numClasses An integer indicating the number of classes to learn
#' @param convergence An float indicating the convergence threshold for stopping iteration
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
    binaryData[fractionIdx]<-rbinom(n=sum(fractionIdx),size=1,prob=data[fractionIdx])
  }

  numSamples=dim(binaryData)[1]            # Number of samples
  posteriorProb=matrix(nrow=numSamples, ncol=numClasses)  # probability each sample (read) belongs to a particular class

  # set.seed(12)    # same random probabilities each time
  for(i in 1:numClasses) {
    # Samples are randomly assigned probabilities (weights) for each class.
    posteriorProb[,i] = rbeta(numSamples,numSamples**-0.5,1)
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
#' @param classMeans
#' @param prev_classMeans
#' @return Returns a vector with the classes ordered
#' @export
order_by_prev_cluster <- function(numClasses, classMeans, prev_classMeans) {
  pr = prev_classMeans
  ord = numeric(length = numClasses)
  for (i in 1:numClasses) {
    # compare cluster means by minimum sum of squares
    pos = which.min(apply(prev_classMeans, 1, function(row) sum((row-classMeans[i,])^2) ))
    ord[i] = pos
    prev_classMeans[pos,] = Inf
  }
  return(ord)
}



#' Classify reads and sort them
#'
#' Classify reads by their posterior probability of belonging to a specific class
#' Then sort the classes by using similarity to the mean profile of previous classes.
#' If the means of previous classes was not provided, hclust is used to cluste classes
#' by their similarity.
#' @param data A matrix of methylation or bincount values (reads x position)
#' @param posteriorProb posteriorProb: a matrix of probabilites of each sample belonging
#'            to a particular class (samples x class)
#' @param previousClassMeans A matrix of the class means from a previous round of clustering
#' @return Returns a matrix with the reads classified (__classX is appended to the read name), and the classes are sorted.
#' @export
classifyAndSortReads<-function(data,posteriorProb,previousClassMeans=NULL) {
  ###################### POST-EM DATA SORTING ###########################################
  #
  # assign classes to reads according to the highest class probability
  numClasses=ncol(posteriorProb)
  readsClasses = apply(posteriorProb, 1, which.max)
  #readsTable = table(readsClasses)
  classMeans = stats::aggregate(data, by = list(readsClasses), FUN = mean)[-1]

  if (!is.null(previousClassMeans)) {
    print("orderByPreviousClusters")
    # order the classes by comparing the class means to the previous clustering
    classOrder = order_by_prev_cluster(numClasses, classMeans, previousClassMeans)
    classMeans = classMeans[classOrder, ]
  } else {
    print("orderByHClust")
    # order the classes by similarity (class means)
    classOrder = stats::hclust(stats::dist(classMeans))$order
    classMeans = classMeans[classOrder, ]
  }

  rownames(data)<-paste(rownames(data),readsClasses,sep="__class")
  ord = order(match(readsClasses,classOrder))
  dataOrderedByClass = data[ord,]

  return(dataOrderedByClass)
}



#' Plot reads by class for a single gene
#'
#' Create a single molecule plot of the reads sorted by class.
#' @param dataOrderedByClass A matrix of methylation or bincount values (reads x position) that have been ordered by class. The assigned class, e.g. "__class1" etc has been appended to read names.
#' @param performSilhouette Make a Silhouette plot to compare within and between class distances
#' @return Produces a single molecule plot sorted by classes and a silouhette plot
#' @export
plotClassesSingleGene<-function(dataOrderedByClass, performSilhouette=T,
                                xlim=c(-250,250), title="Reads by classes",
                                myXlab="CpG/GpC position",
                                featureLabel="TSS", baseFontSize=12) {
  # assign classes to reads according to the highest class probability
  #numClasses=ncol(posteriorProb)
  readClasses <- sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",2)
  readNames<-sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",1)
  readsTable <- table(readsClasses)
  #classMeans = aggregate(data, by = list(readsClasses), FUN = mean)[-1]
  # the horizontal red lines on the plot
  classBorders <- head(cumsum(readsTable[classOrder]), -1)+0.5
  df<-as.data.frame(dataOrderedByClass,stringsAsFactors=F)
  df1<-data.frame(read=readNames,readNumber=1:length(readNames),Class=factor(readClasses),stringsAsFactors=F)

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
                                              col="black")+
      ggplot2::annotate(geom="text", x=1,
                        y=-max(2,0.03*length(reads)),
                        label=featureLabel,color="black")
    # add lines separating classes
    p<-p+ggplot2::geom_hline(yintercept=classBorders,colour="grey80")
    # add color bar for classes
    p<-p+ggplot2::geom_segment(data=df1,mapping=ggplot2::aes(x=xlim[2]+10,y=readNumber-0.5,xend=xlim[2]+10,yend=readNumber+0.5,colour=Class),size=5)+
      ggplot2::geom_vline(xintercept=xlim[2],colour="grey80")
    return(p)
}

  ########################################################################################


#   if (performSilhouette) {
#     region<-sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",1)
#     readClasses <- sapply(strsplit(rownames(dataOrderedByClass),split="__"),"[[",2)
#     df1<-data.frame(regions=region,class=readClasses)
#     table(df1)
#     write.csv(df1,)
#     silhouetteSingleK()
#   }
# }
###############################################################################


