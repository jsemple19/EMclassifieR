% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMmultigene.R
\name{countGenesPerClass}
\alias{countGenesPerClass}
\title{Count genes per class}
\usage{
countGenesPerClass(dataOrderedByClass, sampleName = "")
}
\arguments{
\item{dataOrderedByClass}{A methylation matrix whose row names contain
fields separted by "__" where the first field is the gene name and the
third field is the class number}

\item{sampleName}{String to be used in plot title}
}
\value{
A data frame of counts of reads in each class for each gene
}
\description{
Count number of total genes in multigene clustering and then count
how many unique genes are present for each class (ideally they
should be evenly distributed between the classes)
}
