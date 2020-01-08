# general utility scripts

#' Make directories
#'
#' @param path String with path to where the directories should be made
#' @param dirNameList Vector of strings with names of directories to create (can include multilevel directories)
#' @return Creates the directories listed in dirNameList
#' @examples
#' makeDirs(path=".",dirNameList=c("/txt","/rds/sample1"))
#' @export
makeDirs<-function(path,dirNameList=c()) {
  sub("\\/$","",path) #remove directory slash if present in path string
  for (d in dirNameList) {
    if (!dir.exists(paste0(path,"/",d))){  # for alignments
      dir.create(paste0(path,"/",d), recursive=TRUE, showWarnings=FALSE)
    }
  }
}


#' Get mode (most frequent value)
#'
#' Function copied from https://www.tutorialspoint.com/r/r_mean_median_mode
#' @param v vector of values
#' @return Value that appears most frequently
#' @export
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

