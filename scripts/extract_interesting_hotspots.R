#' -----------------------------------------------------------------------------
#' Extract relevant hotspots based on network properties 
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec  9 09:14:07 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)

dvalidation <- "results/current/biogrid_stringent/validation_tfa/"

# ------------------------------------------------------------------------------
print("meQTL hotspots")
# ------------------------------------------------------------------------------
meqtl <- read_tsv(file.path(dvalidation, "validation_all_meqtl_gwas.txt"))

meqtl_filtered <- meqtl %>% filter(graph_type %in% c("bdgraph", "glasso")) %>%
  filter(!is.na(gwas_disease_trait)) %>%
  filter(snp_cluster_size > 3, snp_genes_selected > 0) %>%
  filter(trans_entities_selected > 0) %>%
  arrange(number_edges, desc(graph_score))

meqtl_top_sentinels <- meqtl_filtered %>% pull(sentinel) %>% unique

# ------------------------------------------------------------------------------
print("eQTL hotspots")
# ------------------------------------------------------------------------------

eqtl <- read_tsv(file.path(dvalidation, "validation_all_eqtlgen_gwas.txt"))

# we use slightly different criteria
eqtl_filtered <- eqtl %>% filter(graph_type %in% c("bdgraph", "glasso")) %>%
  filter(!is.na(gwas_disease_trait)) %>%
  filter(snp_cluster_size > 5, snp_genes_selected > 0, snp_genes_selected < 5) %>%
  filter(trans_entities_selected > 4) %>%
  filter(gwas_trait_match) %>% 
  filter(graph_score > 1) %>%
  arrange(number_edges, desc(graph_score))

# additional filter: one of the selected snp genes matches the best mediating
# gene
eqtl_filtered <- eqtl_filtered %>% 
  separate_rows(snp_genes_selected.list) %>% 
  filter(snp_genes_selected.list == mediation_best_gene)

eqtl_top_sentinels <- eqtl_filtered %>% pull(sentinel) %>% unique

# ------------------------------------------------------------------------------
print("Report results.")
# ------------------------------------------------------------------------------
print("Filtered meQTL sentinels:")
print(meqtl_top_sentinels)
print("Filtered eQTL sentinels:")
print(eqtl_top_sentinels)

print("Sentinels in both datasets:")
print(intersect(meqtl$sentinel, eqtl$sentinel))

# save reduced representaiton of the selected results
meqtl_filtered_reduced <- dplyr::select(meqtl_filtered, sentinel, cohort, graph_type, 
                                        snp_cluster_size, number_edges, graph_score, 
                                        gwas_disease_trait, gwas_trait_match)
eqtl_filtered_reduced <- dplyr::select(eqtl_filtered, sentinel, cohort, graph_type, 
                                       snp_cluster_size, number_edges, graph_score, 
                                       gwas_disease_trait, gwas_trait_match)
write_tsv(meqtl_filtered_reduced,
          file.path(dvalidation, "filtered_meqtl.tsv"))
write_tsv(eqtl_filtered_reduced,
          file.path(dvalidation, "filtered_eqtl.tsv"))

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
