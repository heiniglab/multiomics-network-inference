# ------------------------------------------------------------------------------
#' Script to create the GTEX priors (gene-gene and snp-gene priors)
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# ------------------------------------------------------------------------------
# Load libraries and source scripts
# ------------------------------------------------------------------------------
library(data.table)
library(fdrtool)
source("scripts/priors.R")

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------

# inputs
feqtl <- snakemake@input[["eqtl"]]
dplots <- snakemake@params$plot_dir

# outputs
fout_eqtl_priors <- snakemake@output$eqtl_priors

# ------------------------------------------------------------------------------
print("Start processing.")
# ------------------------------------------------------------------------------
all_priors <- create_eqtlgen_eqtl_priors(feqtl)

# ------------------------------------------------------------------------------
print("Processing done, saving priors.")
# ------------------------------------------------------------------------------
saveRDS(all_priors, fout_eqtl_priors)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
