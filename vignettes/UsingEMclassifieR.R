## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
suppressMessages(library(EMclassifieR))

## -----------------------------------------------------------------------------
tablePath="csv/MatrixLog_relCoord_ampTSS.csv"
matTable<-read.csv(system.file("extdata", tablePath, package = "EMclassifieR",
                               mustWork=TRUE), stringsAsFactors=F)
head(matTable)

## -----------------------------------------------------------------------------
i=1
dataMatrix<-readRDS(system.file("extdata", matTable$filename[i], package = "EMclassifieR", mustWork=TRUE))
dim(dataMatrix)
dataMatrix<-removeNArows(dataMatrix)
dim(dataMatrix)

# dataMatrix<-readRDS("~/Documents/MeisterLab/myPackages/EMclassifieR/vignettes/EMres_cosine/dS02-182_WBGene00015947_K2.rds")
# dataMatrix<-removeNArows(dataMatrix)
# 
# colourChoice=list(low="blue", mid="white", high="red", bg="white", lines="grey80")
# colourChoice=list(low="yellow", mid="grey20", high="red", bg="black", lines="grey80")
# 
# 
# p<-plotClassesSingleMolecule(dataMatrix,
#                                 xRange=c(-250,250), title="Reads by classes",
#                                 myXlab="CpG/GpC position",
#                                 featureLabel="TSS", baseFontSize=12,
#                                 segmentSize=5,
#                                 colourChoice=colourChoice)
# print(p)

## ----singleGeneEM,eval=F------------------------------------------------------
#  
#  k_range = 2:8      # Number of classes to be found
#  maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
#  convergenceError = 10e-6
#  numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
#  xRange=c(-250,250)
#  maxB=50 # Number of randomised matrices to generate
#  outPath="./EMres"
#  maxTasks=4
#  taskId=1
#  nThreads=4
#  setSeed=FALSE
#  distMetric=list(name="cosineDist",rescale=T)
#  
#  if(!dir.exists(outPath)){
#    dir.create(outPath)
#  }
#  
#  #split table indicies into nTasks number of groups
#  taskSubList<-split(1:nrow(matTable),sort(1:nrow(matTable)%%maxTasks))
#  
#  set.seed(200413)
#  for (i in taskSubList[[taskId]]){
#  #for (i in 1:nrow(matTable)) {
#    regionName=matTable$region[i]
#    sampleName=matTable$sample[i]
#    outFileBase=paste(sampleName, regionName, sep="_")
#    print(paste("Clustering", outFileBase))
#    dataMatrix<-readRDS(system.file("extdata", matTable$filename[i],
#                                package = "EMclassifieR", mustWork=TRUE))
#    dim(dataMatrix)
#    dataMatrix<-removeNArows(dataMatrix,maxNAfraction=0.2)
#    #dataMatrix<-recodeMatrixAsNumeric(dataMatrix)
#    dim(dataMatrix)
#  
#    allClassMeans<-tryCatch(
#      {
#        print("running EM for a range of class numbers")
#        runEMrangeClassNum(dataMatrix, k_range, convergenceError, maxIterations,
#                       EMrepeats=numRepeats, outPath=outPath, xRange=xRange,
#                       outFileBase=paste(sampleName, regionName, sep="_"),
#                       doIndividualPlots=FALSE, distMetric=distMetric)
#  
#      },
#        error=function(e){"Matrix not valid"}
#    )
#  
#    if(is.list(allClassMeans)){
#       	saveRDS(allClassMeans,paste0(outPath,"/allClassMeans_",outFileBase,".rds"))
#    } else {
#       	print(allClassMeans) # error message
#    }
#  
#    clustMetrics<-tryCatch(
#      {
#  	    print("plotting clustering metrics for a range of class sizes")
#  	    plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
#  		    maxIterations, outPath, outFileBase, EMrep=NULL, nThreads=nThreads,
#  		    setSeed=setSeed, distMetric=distMetric)
#      },
#      error=function(e){"Matrix not valid"}
#    )
#    if(length(clustMetrics)==1) {
#      print(clustMetrics)
#    }
#  
#    pcaPlots<-tryCatch(
#      {
#        print("plotting PCA of clusters")
#        plotPCAofMatrixClasses(k_range, outPath, outFileBase)
#      },
#      error=function(e){"Matrix not valid"}
#     )
#    if(length(pcaPlots)==1) {
#      print(pcaPlots)
#    }
#  
#    umapPlots<-tryCatch(
#      {
#        print("plotting UMAP of clusters")
#        plotUMAPofMatrixClasses(k_range, outPath, outFileBase)
#      },
#      error=function(e){"Matrix not valid"}
#     )
#    if(length(umapPlots)==1) {
#      print(umapPlots)
#    }
#  
#  }
#  

## ----classesWithTracks,eval=F-------------------------------------------------
#  
#  
#  dataMatrix<-readRDS(system.file("extdata", matTable$filename[i], package = "EMclassifieR", mustWork=TRUE))
#  regionName=matTable$region[i]
#  sampleName=matTable$sample[i]
#  outFileBase=paste(sampleName, regionName, sep="_")
#  
#  allClassMeansList<-readRDS("/Users/semple/Documents/MeisterLab/myPackages/EMclassifieR/inst/extdata/rds/allClassMeans_dS02-182_WBGene00015947.rds")
#  numClasses=5
#  allClassMeans<-allClassMeansList[[numClasses]]
#  allClassMeans
#  
#  ampTSS<-readRDS('/Users/semple/Documents/MeisterLab/myPackages/EMclassifieR/inst/extdata/rds/ampliconMaxTSSgr.RDS')
#  
#  
#  winSize=500
#  tssWin<-ampTSS
#  GenomicRanges::mcols(tssWin)$TSS<-GenomicRanges::start(tssWin)
#  tssWin<-GenomicRanges::resize(tssWin,width=winSize,fix="center")
#  
#  
#  names(GenomicRanges::mcols(tssWin))[1]<-"ID"
#  regionGR<-tssWin[i]
#  featureGR<-ampTSS[i]
#  regionName=matTable$region[i]
#  sampleName=matTable$sample[i]
#  outFileBase=paste(sampleName, regionName, sep="_")
#  #seqlevels(ampTSS)<-seqlevels(Celegans)
#  #seqinfo(ampTSS)<-seqinfo(Celegans)
#  #saveRDS(ampTSS,'/Users/semple/Documents/MeisterLab/myPackages/EMclassifieR/inst/extdata/rds/ampliconMaxTSSgr.RDS')
#  
#  # genomeVer="WS275"
#  # genomeDir=paste0("~/Documents/MeisterLab/GenomeVer/",genomeVer)
#  # txdb<-AnnotationDbi::loadDb(paste0(genomeDir,
#  #                                    "/annotations/c_elegans.PRJNA13758.",
#  #                                    genomeVer,  ".annotations.sqlite"))
#  # seqlevelsStyle(txdb)<-"ucsc"
#  # GenomeInfoDb::genome(txdb)<-GenomeInfoDb::genome(regionGR)
#  #library("TxDb.Celegans.UCSC.ce11.refGene")
#  txdb <- TxDb.Celegans.UCSC.ce11.refGene::TxDb.Celegans.UCSC.ce11.refGene
#  bigwigListFile<-"~/Documents/MeisterLab/sequencingData/public/bigwigFileList.txt"
#  #bigwigListFile<-"~/Documents/MeisterLab/sequencingData/public/bigwigFileList_kranz.txt"
#  strandedBwFile<-"~/Documents/MeisterLab/sequencingData/public/strandedBwFile.txt"
#  gr<-allClassMeansToGR(allClassMeans, regionGR)
#  
#  dTrack<-Gviz::DataTrack(gr,name="dSMF Classes")
#  Gviz::plotTracks(dTrack, groups=rep(c(paste0("class",1:numClasses)), times=10),
#                   type=c("a","p","confint"))
#  
#  print(txdb)
#  outPath="./EMres"
#  grDevices::pdf(paste0(outPath,"/tracks_",
#                            outFileBase,"_K",
#                            numClasses, ".pdf"),
#                     paper="a4", height=11, width=8)
#  
#  plotClassMeansWithTracks(allClassMeans, regionGR, "middle", txdb, featureGR,
#                           bigwigListFile, strandedBwFile)
#  dev.off()

## ----multigeneEM, eval=F------------------------------------------------------
#  tablePath="csv/MatrixLog_relCoord_ampTSS.csv"
#  matTable<-read.csv(system.file("extdata", tablePath, package = "EMclassifieR",
#                                 mustWork=TRUE), stringsAsFactors=F)
#  matTable<-matTable[!is.na(matTable$filename),]
#  
#  head(matTable)
#  multiGeneMat<-NULL
#  genesIncluded<-0
#  for(i in 1:nrow(matTable[matTable$sample=="dS02-182"])){
#    regionName=matTable$region[i]
#    sampleName=matTable$sample[i]
#    outFileBase=paste(sampleName, regionName, sep="_")
#    print(paste0("reading matrix number ",i))
#    dataMatrix<-readRDS(system.file("extdata", matTable$filename[i],
#                                    package = "EMclassifieR", mustWork=TRUE))
#    #Make sure matrix rows do not have too many NAs
#    dataMatrix<-removeNArows(dataMatrix,maxNAfraction=0.2)
#    subMatrix<-selectReadsFromMatrix(dataMatrix,minReads=50,
#                                   addToReadName=outFileBase,
#                                   preferBest=T)
#    if(!is.null(subMatrix)){
#      fullMatrix<-getFullMatrix(subMatrix)
#      winMatrix<-prepareWindows(fullMatrix)
#      genesIncluded<-genesIncluded+1
#      if(is.null(multiGeneMat)){
#        multiGeneMat<-winMatrix
#      } else {
#        multiGeneMat<-rbind(multiGeneMat,winMatrix)
#      }
#    }
#  }
#  print(paste(genesIncluded,"genes included in the multi gene matrix"))
#  
#  #multiGeneMat<-rescale_minus1To1(multiGeneMat)
#  #multiGeneMat<-rescale_0To1(multiGeneMat)
#  #multiGeneMat<-recodeMatrixAsNumeric(multiGeneMat)
#  
#  k_range = 2:8      # Number of classes to be found
#  maxIterations = 100 # number of iterations of EM clustering to perform if it does not converge
#  convergenceError = 10e-6
#  numRepeats=10 # number of repeats of clustering each matrix (to account for fraction of methylation)
#  xRange=c(-250,250)
#  maxB=50 # Number of randomised matrices to generate
#  outPath="./EMres_multigene_cosine_withNAs"
#  maxTasks=4
#  taskId=1
#  nThreads=4
#  setSeed=FALSE
#  distMetric=list(name="cosineDist",rescale=T)
#  
#  if(!dir.exists(outPath)){
#    dir.create(outPath)
#  }
#  
#  
#  set.seed(200413)
#  
#  regionName="multiGene"
#  sampleName="dS02-182"
#  outFileBase=paste(sampleName, regionName, sep="_")
#  print(paste("Clustering", outFileBase))
#  dataMatrix<-multiGeneMat
#  dim(dataMatrix)
#  #dataMatrix<-removeNArows(dataMatrix,maxNAfraction=0)
#  dim(dataMatrix)
#  
#  allClassMeans<-tryCatch(
#    {
#      print("running EM for a range of class numbers")
#      runEMrangeClassNum(dataMatrix, k_range, convergenceError, maxIterations,
#                         EMrepeats=numRepeats, outPath=outPath, xRange=xRange,
#                         outFileBase=paste(sampleName, regionName, sep="_"),
#                         doIndividualPlots=FALSE, distMetric=distMetric)
#  
#    },
#    error=function(e){"Matrix not valid"}
#  )
#  
#  if(is.list(allClassMeans)){
#    saveRDS(allClassMeans,paste0(outPath,"/allClassMeans_",outFileBase,".rds"))
#  } else {
#    print(allClassMeans) # error message
#  }
#  
#  clustMetrics<-tryCatch(
#    {
#      print("plotting clustering metrics for a range of class sizes")
#      plotClusteringMetrics(dataMatrix, k_range, maxB, convergenceError,
#                            maxIterations, outPath, outFileBase, EMrep=NULL,
#                            nThreads=nThreads,
#                            setSeed=setSeed, distMetric=distMetric)
#    },
#    error=function(e){"Matrix not valid"}
#  )
#  if(length(clustMetrics)==1) {
#    print(clustMetrics)
#  }
#  
#  pcaPlots<-tryCatch(
#    {
#      print("plotting PCA of clusters")
#      plotPCAofMatrixClasses(k_range, outPath, outFileBase)
#    },
#    error=function(e){"Matrix not valid"}
#  )
#  if(length(pcaPlots)==1) {
#    print(pcaPlots)
#  }
#  
#  umapPlots<-tryCatch(
#    {
#      print("plotting UMAP of clusters")
#      plotUMAPofMatrixClasses(k_range, outPath, outFileBase)
#    },
#    error=function(e){"Matrix not valid"}
#  )
#  if(length(umapPlots)==1) {
#    print(umapPlots)
#  }
#  
#  print("plotting classes per gene")
#  plotGenesPerClass(k_range, outPath, outFileBase)
#  
#  # to use -1 to 1
#  #https://stats.stackexchange.com/questions/256917/can-an-unstandarized-beta-distribution-have-a-negative-domain
#  

## ----classifyWithKnownClasses,eval=F------------------------------------------
#  allClassMeans<-readRDS(system.file("extdata",
#                                     "EMres/allClassMeans_dS02-182_multiGene.rds",
#                                     package = "EMclassifieR", mustWork=T))
#  repMeans<-allClassMeans[[6]]
#  dd<-repMeans %>% dplyr::group_by(position,class) %>% dplyr::summarize(classMeans=mean(methFreq))
#  dd<-dd %>% tidyr::pivot_wider(names_from=position,values_from=classMeans)
#  dd<-dd[,order(as.numeric(colnames(dd)))]
#  
#  plot(as.vector(dd[1,])~1:481)

