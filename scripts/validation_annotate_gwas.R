#' -----------------------------------------------------------------------------
#' Annotate the validation table with GWAS annotations
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jul 23 13:59:41 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
source("scripts/lib.R")
source("scripts/snipe.R")

fvalidation <- snakemake@input$validation
fgwas <- snakemake@input$gwas

fout_validation <- snakemake@output$validation

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
gwas <- read_tsv(fgwas) %>% as_tibble(.name_repair="universal")

val <- read_tsv(fvalidation)

snps <- unique(val$sentinel)

# get gwas traits for all SNPs
traits_per_snp <- lapply(snps, function(s) {
  return(get_gwas_traits(s, gwas))
})
names(traits_per_snp) <- snps

val$gwas_disease_trait <- NA
val$gwas_disease_trait <- unlist(traits_per_snp[match(val$sentinel,
                                                names(traits_per_snp))])

# ------------------------------------------------------------------------------
print("Writing output file.")
# ------------------------------------------------------------------------------
write_tsv(val, fout_validation, quote_escape=FALSE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

