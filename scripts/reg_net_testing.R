#' -----------------------------------------------------------------------------
#' Script used to unit test methods available in reg_net.R
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

threads <- 4

# ------------------------------------------------------------------------------
print("Load test data.")
# ------------------------------------------------------------------------------
data <- readRDS("results/current/biogrid_stringent/cohort_data_tfa/kora/rs9859077_meqtl.rds")
data <- data.matrix(data)
use <- apply(data,2,function(x) (sum(is.na(x)) / length(x)) < 1)
data <- data[,use]

print(dim(data))

print("testing genenet")
genenet_result <- reg_net(data, NULL, model = "genenet", threads=threads)
print(genenet_result$graph)

print("testing glasso")
glasso_result <- reg_net(data, NULL, model = "glasso", threads = threads)
print(glasso_result$graph)

print("testing genie3")
genie_result <- reg_net(data, NULL, model = "genie3", threads=threads)
print(genie_result$graph)

toplot <- melt(genie_result$fit$pl_fits, id.vars = "weight")
ggplot(toplot, aes(x=weight, y=value)) + facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_line() +
  geom_vline(xintercept = genie_result$fit$best_weight, color="red")


# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
devtools::session_info()

