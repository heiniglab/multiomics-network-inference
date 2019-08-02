print("Installing libraries to:")
.libPaths()

source("https://bioconductor.org/biocLite.R")
biocLite(c("glasso", "GeneNet", "BDgraph", "ROCR", "foreach"))

install.packages("packages/iRafNet_1.1-2.tar.gz", repos=NULL, type="source")
#install.packages("packages/R_peer_source_1.3.tgz", repos=NULL, type="source", INSTALL_opts=c("--no-multiarch"))
