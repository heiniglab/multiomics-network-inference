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

# ------------------------------------------------------------------------------
# Perform validation

# load data and get the validation results
tab <- lapply(ffits, function(f) {
  print(paste0("Processing ", basename(f), "."))

  load(f)

  iteration <- gsub(".RData|.txt","", gsub(".*iter", "", f))
  
  # quick hack to not have to redo output file definitions of other sim parts
  if(any(grepl("prior_completeness", f))) {
    get_validation_table_prior_completeness(result, iteration)
  } else {
    get_validation_table(result, iteration)  
  }
  
}) %>% bind_rows()

# ------------------------------------------------------------------------------
# save results

write_tsv(path=foutput, x=tab)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
