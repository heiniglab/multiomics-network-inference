#' -----------------------------------------------------------------------------
#' Apply the network inference using different models on simulated data and
#' estimate prior completeness
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jan 12 10:45:32 2022
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")

library(GenomicRanges)
library(pheatmap)
library(doParallel)
library(igraph)
library(graph)

source("scripts/lib.R")
source("scripts/reg_net.R")
source("scripts/simulation/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake input and load data.")

# inputs
fdata <- snakemake@input$data
fppi_db <- snakemake@input$ppi_db
fcpg_context <- snakemake@input$cpg_context

# outputs
fout <- snakemake@output[[1]]

# params
threads <- snakemake@threads
sim_iter <- as.numeric(snakemake@params$iteration)

# contains: simulations, ranges, priors, nodes, data, runs
load(fdata)
ppi_db <- readRDS(fppi_db)

# ------------------------------------------------------------------------------
# we generated several graphs, for which we all calculate models now

# apply over the different runs/iterations
run <- simulations[[sim_iter]]
result <- run_ggm_priorcompleteness(run, priors, ranges, 
                                    fcpg_context, ppi_db, threads)

print("Saving results.")
save(file=fout, result)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
