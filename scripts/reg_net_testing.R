#' -----------------------------------------------------------------------------
#' Script used to unit test methods available in reg_net.R
#' Right now no proper unit tests. We merely test whether individual methods
#' throw errors or not.
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Apr  1 12:11:50 2021
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot() + background_grid(major="xy"))

source("scripts/reg_net.R")
source("scripts/lib.R")

threads <- 4

# ------------------------------------------------------------------------------
print("Load test data.")
# ------------------------------------------------------------------------------
data <- readRDS("results/current/biogrid_stringent/cohort_data_tfa/kora/rs9859077_meqtl.rds")
data <- data.matrix(data)
use <- apply(data,2,function(x) (sum(is.na(x)) / length(x)) < 1)
data <- data[,use]

# take a subset of nodes to speed up tests
data <- data[,1:40]
print(dim(data))

priors <- 
  readRDS("results/current/biogrid_stringent/priors/rs9859077_meqtl.rds")
priors <- priors[colnames(data), colnames(data)]

# we set the OMP/BLAS number of threads to 1
# this avoids issues we had in the glasso CV with multi-threading on cluster
# also necessary for BDgraph
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

print("testing genenet")
genenet_result <- reg_net(data, NULL, model = "genenet", threads=threads)
plot(genenet_result$graph)

print("testing irafnet")
irafnet_result <- reg_net(data, priors, model = "irafnet", threads=threads)
plot(irafnet_result$graph)

print("testing glasso")
glasso_result <- reg_net(data, NULL, model = "glasso", threads = threads)
plot(glasso_result$graph)

print("testing glasso with priors")
glasso_result_priors <- reg_net(data, priors, model = "glasso", threads = threads)
plot(glasso_result_priors$graph)

print("testing genie3")
genie_result <- reg_net(data, NULL, model = "genie3", threads=threads)
plot(genie_result$graph)

print("testing bdgraph")
bdgraph_result <- reg_net(data, priors, model = "bdgraph", use_gstart = T,
                          threads=threads)
plot(bdgraph_result$graph)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
devtools::session_info()

