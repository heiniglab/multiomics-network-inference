#' -----------------------------------------------------------------------------
#' Collect probe-ids for all genes available in the eQTL based ranges objects.
#' Also collect rsIDs for the eQTLgen trans-eQTL SNPSs.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Dec 19 20:13:37 2018
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
library(parallel)

# ------------------------------------------------------------------------------
print("Collect all gene probe ids and relevant SNP rsIDs.")
# ------------------------------------------------------------------------------
ppi_db_names <- c("string", "biogrid")
ids <- lapply(ppi_db_names, function(ppi_db_name) {
  franges <-
    list.files(
      paste0("results/current/", ppi_db_name, "/ranges/"),
      "*_eqtlgen.rds$",
      full.names = T
    )
  # probe ids for all genes
  probe_ids <- unique(unlist(lapply(franges, function(fr) {
    ranges <- readRDS(fr)

    c(
      ranges$tfs$ids,
      ranges$spath$ids,
      ranges$trans_genes$ids,
      ranges$snp_genes$ids
    )

  })))

  # rs ids for all trans eQTL SNPs
  rs_ids <- unique(unlist(lapply(franges, function(fr) {
    ranges <- readRDS(fr)

    names(ranges$sentinel)
  })))
  write.table(file=paste0("eqtlgen_probe_ids_", ppi_db_name, ".txt"),
              cbind(probe_ids=probe_ids), col.names=F, row.names=F, quote=F)
  write.table(file=paste0("eqtlgen_rs_ids_", ppi_db_name, ".txt"),
              cbind(rs_id=rs_ids), col.names=F, row.names=F, quote=F)

  list(probe_ids=probe_ids, rs_ids=rs_ids)
})
names(ids) <- ppi_db_names

# compare with available lolipop data
load("results/current/ggmdata_lolipop.RData")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
