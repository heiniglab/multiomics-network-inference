#' -----------------------------------------------------------------------------
#' Summarize SNP gene recovery in simulation for a specific sentinel
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed May  6 11:23:47 2020
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
library(graph)

# the sentinel to look for
sentinel <- snakemake@wildcards$sentinel

# the directory containing the simulation results
dresults <- snakemake@params$dresults

# the output file
fout <- snakemake@output$summary

threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Load and process data.")
# ------------------------------------------------------------------------------
library(parallel)

# we use all results from the "full" sample size simulation
finputs <- list.files(dresults, paste0(sentinel, ".*-subsetall.RData"), 
                      full.names = T)

if(length(finputs) < 100) stop("Missing results.")

res <- mclapply(finputs, function(f) {
  print(paste0("Processing ", f, "."))
  
  n <- load(f) # result object
  if(length(result) == 0) {warning("no results"); return(NULL)}
  
  # only look at results with full prior information
  result <- result[grepl("_rd0$", names(result))][[1]]
  
  # get observed graph and extract SNP genes
  gobs <- result$graph.observed
  snp <- result$snp
  if(!snp %in% nodes(gobs)) return(NULL)
  
  sgenes <- adj(gobs, result$snp)[[1]]
  
  if(length(sgenes) == 0) return(NULL)
  
  # iteratie over all inferred graphs and get number of correctly inferred 
  # snp genes
  inferred_graphs <- result$fits[!grepl("_fit$", names(result$fits))]
  
  lapply(names(inferred_graphs), function(n) {
    ginf <- inferred_graphs[[n]]
    
    if(!snp %in% nodes(ginf)) return(NULL)
    
    sgenes_inf <- adj(ginf, result$snp)[[1]]
    
    # get stats
    tp <- length(intersect(sgenes_inf, sgenes))
    fp <- length(setdiff(sgenes_inf, sgenes))
    acc <- tp / length(sgenes)
    
    tibble(result$snp, graph_type = n, number_snp_genes = length(sgenes),
           recovered_snp_genes = length(sgenes_inf),
           TP = tp, FP = fp, ACC = acc)
    
  }) %>% bind_rows()
}, mc.cores = threads) %>% bind_rows()

# plot the results and save table
write_tsv(res, fout)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
