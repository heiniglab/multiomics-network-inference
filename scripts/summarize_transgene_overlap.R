#' -----------------------------------------------------------------------------
#' description
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Jan 18 21:15:49 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(ggplot2)

# ------------------------------------------------------------------------------
print("Loading and preprocessing data.")
# ------------------------------------------------------------------------------
#output
fout <- snakemake@output[[1]]

# params
hotspot_threshold <- as.numeric(snakemake@params$hotspot_threshold)

feqtl <- snakemake@input[[1]]
eqtl <- readRDS(feqtl)
eqtl$SNPChr <- paste0("chr", eqtl$SNPChr)

# split by chromosome
eqtl_by_chromosome <- split(eqtl, f=eqtl$SNPChr)

overlaps <- lapply(names(eqtl_by_chromosome), function(chr) {
  eqtl_subs <- eqtl_by_chromosome[[chr]]
  eqtl_subs <- eqtl_subs[order(eqtl_subs$SNPPos),]

  # index for hotspots
  filter <- tapply(eqtl_subs$Gene, eqtl_subs$SNP, length) > hotspot_threshold
  filter <- filter[filter]
  hotspots <- eqtl_subs[eqtl_subs$SNP %in% names(filter),]

  # same for pruned
  pruned_subs <- eqtl_subs[eqtl_subs$SNP == eqtl_subs$sentinel,]
  filter <- tapply(pruned_subs$Gene, pruned_subs$SNP, length) > hotspot_threshold
  filter <- filter[filter]
  hotspots_pruned <- eqtl_subs[eqtl_subs$SNP %in% names(filter),]

  get_gene_overlaps <- function(hs) {
    sentinels <- unique(hs$SNP)

    overlaps <- list()
    for(i in 1:length(sentinels)) {
      if(i == length(sentinels)) next;

      s1 <- sentinels[i]
      s2 <- sentinels[i+1]
      genes1 <- hs[hs$SNP == s1, "GeneSymbol"]
      genes2 <- hs[hs$SNP == s2, "GeneSymbol"]

      overlap <- intersect(genes1, genes2)
      overlap_ratio <- length(overlap) / max(length(genes1), length(genes2))

      overlaps[[i]] <- c(length(overlap), overlap_ratio)
    }
    return(overlaps)
  }
  go <- get_gene_overlaps(hotspots)
  go_pruned <- get_gene_overlaps(hotspots_pruned)

  list(full=matrix(unlist(go), ncol=2, byrow = T),
       pruned=matrix(unlist(go_pruned), ncol=2, byrow = T))
})

# get individual matrices for pruned and full sets
overlaps_pruned <- do.call(rbind, lapply(overlaps, "[[", "pruned"))
overlaps_full <- do.call(rbind, lapply(overlaps, "[[", "full"))

colnames(overlaps_pruned) <- colnames(overlaps_full) <- c("overlap_count",
                                                          "overlap_ratio")

# create data.frame for plotting
toplot <- rbind(overlaps_full,overlaps_pruned)
toplot <- cbind.data.frame(toplot, group = c(rep("full", nrow(overlaps_full)),
                                  rep("pruned", nrow(overlaps_pruned))))

# plot distributions
pdf(fout)
theme_set(theme_bw())
ggplot(toplot, aes(x=overlap_count, color=group)) +
  geom_freqpoly() +
  ggtitle("Distribution of number of trans-gene overlaps for neighbouring hotspots\nbefore and after SNP pruning")
ggplot(toplot, aes(x=overlap_ratio, color=group)) +
  geom_freqpoly() +
  ggtitle("Distribution of fraction of trans-gene overlap for neighbouring hotspots\nbefore and after SNP pruning")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
