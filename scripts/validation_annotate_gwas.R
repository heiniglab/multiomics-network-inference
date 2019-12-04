#' -----------------------------------------------------------------------------
#' Annotate the validation table with GWAS annotations
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jul 23 13:59:41 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

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
  as_tibble(.name_repair="universal") %>%
  separate_rows(SNPS)

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
  eqtl %>% 
    filter(FDR < 0.01 & GeneSymbol %in% gene_list) %>%
    pull(SNP) %>% unique
}
print("Loading eQTLgen data.")
eqtlgen <- read_tsv(feqtlgen)

print("Procesing without SNiPA results.")
traits_per_locus_genes <- mclapply(val$non_snp_genes_selected.list, function(genes) {
  if(!is.na(genes)) {
    gene_list <- strsplit(genes, ";")[[1]]
    
    snps <- sapply(gene_list, function(g) {
      s <- get_eqtl_snps(g, eqtlgen)  
      
      # only keep rs-ids, otherwise SNiPA will not work
      s[grepl("^rs", s)]
    })
    
    traits <- get_gwas_traits(unique(unlist(snps)), gwas, 
                              get.ld.snps = F, collapse = F)
    paste0(unique(traits), collapse = ";")
    
  } else {
    NA
  }
}, mc.cores=threads)
val$non_snp_genes_selected.gwas_traits <- unlist(traits_per_locus_genes)
print("Done.")

# any of the SNP traits matches one of the gene traits?
print("Checking SNP/Gene trait matches.")
val$gwas_trait_match <- sapply(1:nrow(val), function(i) {
  straits <- val$gwas_disease_trait[i]
  gtraits <- val$non_snp_genes_selected.gwas_traits[i]
  
  out <- FALSE
  if(!is.na(straits) & !is.na(gtraits)) {
    straits <- strsplit(straits, "\\|")[[1]]
    # split up gene traits, too, to make sure we match complete trait names
    gtraits <- strsplit(gtraits, ";")[[1]]
    if(any(sapply(straits, function(tr) { any(grepl(tr, gtraits)) }))) {
      out <- TRUE
    }
  }
  out
})

# ------------------------------------------------------------------------------
print("Writing output file.")
# ------------------------------------------------------------------------------
write_tsv(val, fout_validation)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

