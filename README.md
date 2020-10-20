# EMclassifieR
Package for classifying methylation matrices from reads into probabilistic classes using the expectation maximisation algorithm. Based on work by Kyriakos Schwartz, which was based on the following publication, and code in the supplementary material thereof:

__Nair, N. U., Kumar, S., Moret, B. M. E., & Bucher, P. (2014). Probabilistic partitioning methods to find significant patterns in ChIP-Seq data. Bioinformatics, 30(17), 2406–2413. https://doi.org/10.1093/bioinformatics/btu318__

Kyriakos Schwartz did an initial implementation of the code for binary methylation data. Here I turn this code in to an R package and expand its functionality.

# Installation
Install some necessary packages:
```
install.packages(c("zoo", "cluster", "ggplot2", "tidyr", "plyr", "purrr", "ggpubr", "parallel", "foreach", "doParallel", "magrittr", "doRNG", "lsa", "ggbiplot", "umap", "HiClimR"))

BiocManager::install(c("GenomicRanges","IRanges","GenomeInfoDb","rtracklayer","Gviz","TxDb.Celegans.UCSC.ce11.refGene","maigesPack"))
devtools::install_github("vqv/ggbiplot")
```

Then install this package:
```
devtools::install_github("jsemple19/EMclassifieR")
```
