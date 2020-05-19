#' -----------------------------------------------------------------------------
#' Extract SNP info from the Lean Body Mass GWAS
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Feb  5 07:13:10 2020
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
source("scripts/biomaRt.R")

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------
fgwas <- snakemake@input$gwas
fout_snp_locs <- snakemake@output$snp_locs
fout_snp_pvalues <- snakemake@output$snp_pvalues

# load GWAS data and Position data for each SNP in there
gwas <- read_tsv(fgwas)
snp_pos <- get_snpPos_biomart(gwas$MarkerName)

# merge GWAS with position data, prepare dataframe for BIM file output
gwas_with_pos <- left_join(gwas, snp_pos, by=c("MarkerName" = "snp")) %>% 
  filter(!is.na(chr)) %>% 
  mutate(filler=0, A1=toupper(Allele1), A2=toupper(Allele2)) %>% 
  dplyr::select(chr, snp=MarkerName, filler, start, A1, A2, `P-value`) %>% 
  filter(!grepl("^H", chr))

# write bim output format
write_tsv(gwas_with_pos %>% dplyr::select(-`P-value`), 
          fout_snp_locs,
          col_names=F)

# write P-value file
write_tsv(gwas_with_pos %>% dplyr::select(snp, `P-value`), 
          fout_snp_pvalues,
          col_names=T)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
