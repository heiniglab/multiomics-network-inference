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
library(parallel)
source("scripts/lib.R")
source("scripts/snipe.R")

fvalidation <- snakemake@input$validation
fgwas <- snakemake@input$gwas
feqtlgen <- snakemake@input$eqtlgen

fout_validation <- snakemake@output$validation

threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
gwas <- read_tsv(fgwas, col_types = cols(.default="c")) %>%
  as_tibble(.name_repair="universal")

# avoid the cluster_sizes column to be parsed as an integer (it's a comma
# separated string)
val <- read_tsv(fvalidation, col_types = cols(.default="c"))

snps <- unique(val$sentinel)

# get gwas traits for all SNPs
traits_per_snp <- mclapply(snps, function(s) {
  return(get_gwas_traits(s, gwas))
}, mc.cores=threads)
names(traits_per_snp) <- snps

val$gwas_disease_trait <- NA
val$gwas_disease_trait <- unlist(traits_per_snp[match(val$sentinel,
                                                names(traits_per_snp))])

# also process the gene list: get cis-eQTL snps and check their GWAS annot
print("Annotating selected genes with eQTL SNP GWAS traits.")
get_eqtl_snps <- function(gene_list, eqtl) {
  eqtl %>% filter(GeneSymbol %in% gene_list) %>%
    pull(SNP) %>% unique
}
print("Loading eQTLgen data.")
eqtlgen <- read_tsv(feqtlgen)

print("Annotating...")
traits_per_locus_genes <- mclapply(val$non_snp_genes_selected.list, function(genes) {
  if(!is.na(genes)) {
    gene_list <- strsplit(genes, ";")[[1]]
    results <- sapply(gene_list, function(g) {
      snps <- get_eqtl_snps(g, eqtlgen)  
      
      # only keep rs-ids, otherwise SNiPA will not work
      snps <- snps[grepl("^rs", snps)]
      get_gwas_traits(snps, gwas)
    })
    paste0(gene_list, ":", results, collapse = ";")
    
  } else {
    NA
  }
}, mc.cores=threads)
print("Done.")
val$non_snp_genes_selected.gwas_traits <- unlist(traits_per_locus_genes)

sapply(val$gwas)
# ------------------------------------------------------------------------------
print("Writing output file.")
# ------------------------------------------------------------------------------
write_tsv(val, fout_validation)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

