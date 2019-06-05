print("Installing libraries to:")
.libPaths()

source("https://bioconductor.org/biocLite.R")
biocLite(c("glasso", "GeneNet", "BDgraph","igraph", "ROCR", "doParallel", "foreach"))

install.packages("packages/iRafNet_1.1-2.tar.gz", repos=NULL, type="source")
