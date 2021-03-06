% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessMethMats.R
\name{removeNArows}
\alias{removeNArows}
\title{Remove rows containing NAs}
\usage{
removeNArows(dataMatrix, maxNAfraction = 0.2, removeAll0 = F)
}
\arguments{
\item{dataMatrix}{A matrix of numbers and NAs.}

\item{maxNAfraction}{Maximual fraction of CpG/GpC positions that can be undefined (default=0.2)}

\item{removeAll0}{Remove reads that only have 0 or NA values, i.e. were not methylated at all (default=F)}
}
\value{
A matrix without the rows where the fraction of NA positions is above
the threshold
}
\description{
Function to remove all rows from a matrix where more than a certain
fraction of the positions in that row are NAs. To remove rows with any
NAs at all, set maxNAfraction=0
}
\examples{
removeNArows(matrix(c(1,2,3,NA,4,5),nrow=3,byrow=TRUE))
}
