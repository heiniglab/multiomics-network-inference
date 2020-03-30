# ------------------------------------------------------------------------------
#'
#' Script to validate the GGMs calculated on simulated data
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
sink(file=snakemake@log[[1]])

# ------------------------------------------------------------------------------
# load needed libraries and source scripts

library(BDgraph)
library(graph)
library(igraph)
library(GenomicRanges)
library(tidyverse)
source("scripts/lib.R")
source("scripts/simulation/lib.R")

# ------------------------------------------------------------------------------
# get snakemake params

# inputs
ffits <- snakemake@input$fits

# outputs
foutput <- snakemake@output[[1]]

# iteration/run to be annotated in result table
iteration <- snakemake@params$iteration

print(paste0("Processing file: ", ffits, "."))

# ------------------------------------------------------------------------------
# Perform validation

# load data and get the validation results
load(ffits)
tab <- get_validation_table(result, iteration)

# ------------------------------------------------------------------------------
# save results

write_tsv(path=foutput, x=tab)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
