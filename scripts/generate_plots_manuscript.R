#' -----------------------------------------------------------------------------
#' Script to generate the individual plots we use in the manuscript.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Oct 21 18:22:43 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(grid)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
library(circlize)
library(ggbio)
library(xtable)

# get the hg19 chromosome definitions
library(BSgenome.Hsapiens.UCSC.hg19)
hg19info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

# set up theme and colors
theme_set(theme_cowplot() + background_grid())
theme_update(legend.text = element_text(size=11), 
             legend.title=element_text(size=12),
             axis.text.x = element_text(size=10),
             axis.text.y = element_text(size=10))

# manually define the colors for the different methods (to get paired ones for
# prior/non prior versions)
paired <- brewer.pal(4, "Paired")
names(paired) <- c("glasso", "glasso (priors)", "bdgraph", "bdgraph (priors)")
#paired <- c(paired, "#66D35F")
#names(paired) <- c("glasso", "glasso (priors)", 
#                   "bdgraph (empty)", "bdgraph (priors)", "bdgraph (full)")
unpaired <- brewer.pal(7, "Dark2")[c(2,3,7)]
names(unpaired) <- c("irafnet", "genie3", "genenet")
graph_cols <- c(paired, unpaired)

sfb_graphs <- scale_fill_manual(values=graph_cols)
#sfb_graphs <- scale_fill_brewer(palette="Set2")
#scb_graphs <- scale_color_brewer(palette="Set2")
scb_graphs <- scale_color_manual(values=graph_cols)
sfb_binary <- scale_fill_brewer(palette = "Accent")
scb_binary <- scale_color_brewer(palette = "Accent")
scb_priors <- scale_color_brewer(palette = "Dark2")

bgm <- background_grid(major = "xy")
group_cols <- brewer.pal("Set2", n=3)
COLORS <- list(MEQTL = group_cols[1],
               EQTL = group_cols[2])

# for the main results (not including hotspot information)
RESULT_PATH <- "results/current/biogrid_stringent/"

# ------------------------------------------------------------------------------
print("Figure 1 - Panel A")
# ------------------------------------------------------------------------------
# This panel is the graphical abstract of the study and was manually designed 
# in inkscape.

# ------------------------------------------------------------------------------
print("Figure 1 - Panel B")
# ------------------------------------------------------------------------------
eqtl_hotspots <- read_tsv("results/current/hotspots/eqtlgen_thres5/hotspots.tsv")

eqtl_regions <- with(eqtl_hotspots, 
                     GRanges(paste0("chr", SNPChr),
                             IRanges(SNPPos, width=1),
                             name=SNP,
                             trans_associations=ntrans,
                             seqinfo = hg19info))

eqtl_regions$trans_ranges <- with(eqtl_hotspots, 
                                  GRanges(paste0("chr", GeneChr),
                                          IRanges(GenePos, width=2),
                                          name=GeneSymbol,
                                          seqinfo=hg19info))

meqtl_hotspots <- read_tsv("results/current/hotspots/meqtl_thres5/hotspots.tsv")

meqtl_regions <- with(meqtl_hotspots, 
                      GRanges(paste0("chr", chr.snp),
                              IRanges(interval.start.snp + (interval.end.snp-interval.start.snp)/2, 
                                      width=1),
                              name=sentinel.snp,
                              trans_associations=ntrans,
                              seqinfo = hg19info))

meqtl_regions$trans_ranges <- with(meqtl_hotspots, 
                                   GRanges(paste0("chr", chr.cpg),
                                           IRanges(interval.start.cpg + (interval.end.cpg-interval.start.cpg)/2, width=2),
                                           name=sentinel.cpg))


# Plot QTL pairs. Code adapted from Julian Schmidt
plot_pairs <- function(genome, ranges_x, ranges_y, log=FALSE, 
                       label_x="x", label_y="y", resolution=1e7,
                       color = brewer.pal("Reds", n=3)[1], ...) {
  chrs = paste0("chr", 1:22)
  chrlen <- seqlengths(genome)
  chrlen <- chrlen[chrs]
  
  genome_bins <- tileGenome(chrlen, tilewidth = resolution, cut.last.tile.in.chrom = T)
  
  breaks <- table(seqnames(genome_bins))
  
  for (i in 2:length(breaks)) {
    breaks[i] <- breaks[i - 1] + breaks[i]
  }
  
  bin_overlaps_x <- findOverlaps(ranges_x, genome_bins)
  bin_overlaps_y <- findOverlaps(ranges_y, genome_bins)
  
  mappable_pairs <-
    intersect(queryHits(bin_overlaps_x), queryHits(bin_overlaps_y))
  
  x_bin <-
    subjectHits(bin_overlaps_x)[queryHits(bin_overlaps_x) %in% mappable_pairs]
  y_bin <-
    subjectHits(bin_overlaps_y)[queryHits(bin_overlaps_y) %in% mappable_pairs]
  
  pairs_binned <- cbind.data.frame(x_bin, y_bin)
  pairs_binned <-
    pairs_binned[order(pairs_binned$x_bin, pairs_binned$y_bin), ]
  
  x_margin <- group_by(pairs_binned, x_bin) %>% summarise(count=n())
  y_margin <- group_by(pairs_binned, y_bin) %>% summarise(count=n())
  
  pairs_binned_discrete <- pairs_binned
  
  pairs_binned <- data.frame(table(pairs_binned), stringsAsFactors=F)
  pairs_binned[,1] <- as.numeric(as.character(pairs_binned[,1]))
  pairs_binned[,2] <- as.numeric(as.character(pairs_binned[,2]))
  pairs_binned[pairs_binned$Freq > 50, "Freq"] <- 50
  
  xmp <- ggplot(x_margin) + 
    geom_density(aes(x=x_bin, y=count), stat="identity", color=color, shape=23) + 
    scale_x_continuous(expand = c(0.01, 0.01), breaks = as.vector(breaks), labels = NULL) + 
    xlab("") + 
    background_grid(major = "xy") + 
    theme(axis.text.y = element_text(size=9),
          axis.title.y = element_text(size=10))
  
  ymp <- ggplot(y_margin) + 
    geom_density(aes(x=y_bin, y=count), stat="identity", color=color, shape=23) + 
    coord_flip() + scale_x_continuous(
      expand = c(0.01, 0.01),
      breaks = as.vector(breaks),
      labels = NULL) + 
    xlab("") + 
    theme(axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=9,
                                     angle=90, 
                                     vjust=0.5, 
                                     hjust=1)) + 
    background_grid(major = "xy")
  
  # the discretized plot
  gd <-
    ggplot(pairs_binned_discrete) + 
    geom_tile(..., fill=color,
              aes(x = x_bin, y = y_bin)) +
    theme(
      text = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(size=9, vjust = 0.5, angle = 90),
      axis.text.y = element_text(size=9),
      legend.title = element_text(size = 10),
      plot.margin = margin(0,0.1,0,0, "cm")) +
    xlab(label_x) + ylab(label_y) +
    scale_x_continuous(
      expand = c(0.01, 0.01),
      breaks = as.vector(breaks),
      labels = names(breaks),
      limits = c(1,length(genome_bins))
    ) +
    scale_y_continuous(
      expand = c(0.01, 0.01),
      breaks = as.vector(breaks),
      labels = names(breaks)
    ) + background_grid(major = "xy")
  
  list(x_margin_plot=xmp, y_margin_plot=ymp, data=pairs_binned, 
       data_disc = pairs_binned_discrete, pair_plot_disc=gd)
}

tile_width <- 5

eqtl_plot <- plot_pairs(hg19info, 
                        eqtl_regions, 
                        eqtl_regions$trans_ranges, 
                        log = F, "SNPs", "trans Genes", resolution = 2e6, width=tile_width, color=COLORS$EQTL)
meqtl_plot <- plot_pairs(hg19info, 
                         meqtl_regions, 
                         meqtl_regions$trans_ranges, 
                         log=F, "SNPs", "CpGs", resolution = 5e6, width=tile_width, color=COLORS$MEQTL)

wr <- c(85,15)
hr <- c(15,85)
side_margins <- margin(0.1,0.7,0.1,-0.9, unit = "lines")
top_margins <- margin(0.7, .3,-0.8, 0.5, unit = "lines")

# plot the discretized version
# meQTL subplot
ap <- meqtl_plot$x_margin_plot
bp <- meqtl_plot$pair_plot_disc
cp <- meqtl_plot$y_margin_plot
ap <- ap + theme(plot.margin = top_margins)
cp <- cp + theme(plot.margin = side_margins)

# combine our individual grobs
meqtl_ggplot <- ggarrange(ggarrange(ap, nullGrob(), widths = wr, ncol=2), 
                          ggarrange(bp, cp, widths=wr, ncol=2, align="h"), 
                          nrow=2, heights = hr, align="v")

# eqtl subplot
top_margins <- margin(0.7, 0.2,-0.8, 0.1, unit = "lines")
ap <- eqtl_plot$x_margin_plot
bp <- eqtl_plot$pair_plot_disc
cp <- eqtl_plot$y_margin_plot
ap <- ap + theme(plot.margin = top_margins)
cp <- cp + theme(plot.margin = side_margins)

# combine again
eqtl_ggplot <- ggarrange(ggarrange(ap, nullGrob(), widths = wr, ncol=2), 
                         ggarrange(bp, cp, widths=wr, ncol=2, align="h"), 
                         nrow=2, heights = hr, align="v")

panel_b1 <- meqtl_ggplot
panel_b2 <- eqtl_ggplot

# used in the thesis
ga <- ggarrange(panel_b1, panel_b2, ncol=2, labels = "AUTO")
save_plot("results/current/figures/thesis/hotspots_pairplot.pdf", ga, ncol=3, nrow=2)

# ------------------------------------------------------------------------------
print("Figure 1 - Panel B alternative using grandLinear plots, currently used")
# ------------------------------------------------------------------------------
hg19info2 <- dropSeqlevels(keepStandardChromosomes(hg19info), c("chrX", "chrY", "chrM"))

eqtl_hotspots <- read_tsv("results/current/hotspots/eqtlgen_thres5/hotspots.tsv")

eqtl_regions <- unique(with(eqtl_hotspots, 
                     GRanges(paste0("chr", SNPChr),
                             IRanges(SNPPos, width=1),
                             name=SNP,
                             trans_associations=ntrans,
                             seqinfo = hg19info2)))

meqtl_hotspots <- read_tsv("results/current/hotspots/meqtl_thres5/hotspots.tsv") %>%
  arrange(chr.snp) %>%
  filter(sentinel.snp == snp)

meqtl_regions <- unique(with(meqtl_hotspots, 
                      GRanges(paste0("chr", chr.snp),
                              IRanges(pos.snp, 
                                      width=1),
                              name=sentinel.snp,
                              trans_associations=ntrans,
                              seqinfo = hg19info2)))

eqtl_regions$category <- "eQTL"
meqtl_regions$category <- "meQTL"
regions <- c(eqtl_regions, meqtl_regions)

gl_eqtl <- plotGrandLinear(eqtl_regions, aes(y=trans_associations), shape=2, size=3,
                color=COLORS$EQTL, coord="genome", spaceline = T) + 
  theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0, size=9)) + 
  labs(y = "Number of trans-eGenes")

gl_meqtl <- plotGrandLinear(meqtl_regions, aes(y=trans_associations), shape=2, size=3,
                color=COLORS$MEQTL, coord="genome", spaceline = T) + 
  theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0, size=9)) + 
  labs(y = "Number of trans-CpGs")

panel_b1 <- gl_eqtl@ggplot 
panel_b2 <- gl_meqtl@ggplot

# ------------------------------------------------------------------------------
print("Figure 1 - Panel C")
# ------------------------------------------------------------------------------
# get meqtl ranges
fmeqtl_ranges <- list.files(paste0(RESULT_PATH, "ranges/"), "*_meqtl.rds", 
                            full.names = T)

# extract entity counts
meqtl_counts <- lapply(fmeqtl_ranges, function(fi) {
  ranges <- readRDS(fi)
  
  # filter out genes for which we didnt have and probe ids
  # (those were not used in inference)
  ncpgs <- length(ranges$cpgs)
  ntfs <- length(ranges$tfs[!sapply(ranges$tfs$ids, is.null)])
  nspath <- length(ranges$spath[!sapply(ranges$spath$ids, is.null)])
  nsnp_genes <- length(ranges$snp_genes[!sapply(ranges$snp_genes$ids, is.null)])
  ncpg_genes <- length(ranges$cpg_genes[!sapply(ranges$cpg_genes$ids, is.null)])
  
  data.frame(cis_genes = nsnp_genes, 
             trans_associations = ncpgs,
             trans_genes = ncpg_genes,
             TFs = ntfs,
             PPI = nspath, 
             group = "meQTL")
})
meqtl_counts <- bind_rows(meqtl_counts)

# get the eqtl entity counts
feqtl_ranges <- list.files(paste0(RESULT_PATH, "ranges/"), "*_eqtlgen.rds", 
                           full.names = T)
eqtl_counts <- lapply(feqtl_ranges, function(fi) {
  ranges <- readRDS(fi)
  
  # filter out genes for which we didnt have and probe ids
  # (those were not used in inference)
  ntrans_genes <- length(ranges$trans_genes[!sapply(ranges$trans_genes$ids, is.null)])
  ntfs <- length(ranges$tfs[!sapply(ranges$tfs$ids, is.null)])
  nspath <- length(ranges$spath[!sapply(ranges$spath$ids, is.null)])
  nsnp_genes <- length(ranges$snp_genes[!sapply(ranges$snp_genes$ids, is.null)])
  
  data.frame(cis_genes = nsnp_genes, 
             trans_associations = ntrans_genes,
             trans_genes = NA,
             TFs = ntfs,
             PPI = nspath,
             group = "eQTL")
})
eqtl_counts <- bind_rows(eqtl_counts)

# gather in single df
counts <- bind_rows(meqtl_counts, eqtl_counts)

# get some numbers for the manuscript, too
medians <- function(d, g) {
  d %>%
    filter(group %in% g) %>% 
    dplyr::select(-group) %>%
    summarise_all(.funs = median)
}
print("Median number of entities for eQTLs:")
medians(counts, "eQTL")
print("Median number of entities for meQTLs:")
medians(counts, "meQTL")

# create the actual plot
toplot <- melt(counts) %>% 
  arrange(value) %>% 
  mutate(variable = gsub("cis_genes", "cis genes", variable)) %>%
  mutate(variable = gsub("trans_associations", "trans CpGs/genes", variable)) %>%
  mutate(variable = gsub("trans_genes", "CpG genes", variable)) %>%
  mutate(variable = gsub("PPI", "PPI genes", variable))

panel_c <- ggplot(toplot, aes(color=group, x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=23, position = position_jitterdodge(jitter.width = 0.15)) + 
  scale_color_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL)) + 
  scale_y_log10() +
  labs(x="",
       y="count") + 
  theme(legend.position = "none", 
        plot.margin = margin(1,0.5, 0, 0.2, unit="lines"),
        axis.text.x = element_text(size=11,angle=-45, hjust=0,vjust=1))


# generate the plot for the thesis
panel_c1_thesis <- toplot %>%
  filter(group == "meQTL") %>%
  ggplot(aes(color=group, x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=23, position = position_jitterdodge(jitter.width = 0.15)) + 
  scale_color_manual(values=c(COLORS$MEQTL)) + 
  scale_y_log10() +
  labs(x="",
       y="count") + 
  theme(legend.position = "none", 
        plot.margin = margin(1,0.5, 0, 0.2, unit="lines"),
        axis.text.x = element_text(size=11,angle=-45, hjust=0,vjust=1))

panel_c2_thesis <- toplot %>%
  filter(group == "eQTL") %>%
  drop_na(value) %>%
  ggplot(aes(color=group, x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=23, position = position_jitterdodge(jitter.width = 0.15)) + 
  scale_color_manual(values=c(COLORS$EQTL)) + 
  scale_y_log10() +
  labs(x="",
       y="count") + 
  theme(legend.position = "none", 
        plot.margin = margin(1,0.5, 0, 0.2, unit="lines"),
        axis.text.x = element_text(size=11,angle=-45, hjust=0,vjust=1))

panel_c_thesis <- ggarrange(panel_c1_thesis,
                            panel_c2_thesis,
                            ncol=2, labels="AUTO")

save_plot("results/current/figures/thesis/collected_entities_thesis.pdf", panel_c_thesis,
          ncol=2, base_asp = 1.2)

# ------------------------------------------------------------------------------
print("Figure 1 - Panel D")
# ------------------------------------------------------------------------------
# prior plot
finput <- list.files(paste0(RESULT_PATH, "priors/"), "*.rds", full.names = T)
tab <- lapply(finput, function(f) {
  priors <- readRDS(f)

  total_priors <- unname(table(priors[upper.tri(priors)]>min(priors))["TRUE"])
  total_nodes <- ncol(priors)
  total_edges <- (total_nodes * (total_nodes-1)) / 2
  if(grepl("meqtl", f)) {
    group = "meQTL"
  } else {
    group = "eQTL"
  }
  c(total_priors, total_nodes, total_edges, group)
})

tab <- cbind.data.frame(total_priors = as.numeric(sapply(tab, "[[", 1)),
                        total_nodes = as.numeric(sapply(tab, "[[", 2)),
                        total_edges = as.numeric(sapply(tab, "[[", 3)),
                        group = sapply(tab, "[[", 4),
                        stringsAsFactors=F)
tab$fraction_priors <- tab$total_priors / tab$total_edges

# this might be an interesting option for the paper.
# Another ideas is to stratify the priors by prior type and create a freqpoly plot
panel_d <- tab %>%
  ggplot(aes(x=total_priors)) + 

  geom_freqpoly(aes(color=group), size=1, bins = 100) +
  scale_color_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL)) + 
  geom_vline(xintercept = min(tab$total_priors), col="#777777", linetype="dashed") + 
  annotate(geom="text", 
           x=min(tab$total_priors),
           y=0, 
           size=2.5,
           label=paste0("min(x) = ", min(tab$total_priors)), 
           hjust=-0.1, 
           vjust=1.2, 
           col="#777777") + 
  coord_cartesian(xlim=c(0,1000)) +
  labs(x="number of edges with priors") + 
  theme(plot.margin = margin(1,1.5,1,1, unit="lines"), 
        legend.position = c(0.7,0.9),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

# some stats
median(tab %>% as_tibble() %>% filter(group=="meQTL") %>% dplyr::pull(total_priors))
median(tab %>% as_tibble() %>% filter(group=="eQTL") %>% dplyr::pull(total_priors))
length(tab %>% as_tibble() %>% filter(group=="eQTL" & total_priors >= 100) %>% arrange(total_priors) %>% dplyr::pull(total_priors))
length(tab %>% as_tibble() %>% filter(group=="meQTL" & total_priors >= 100) %>% arrange(total_priors) %>% dplyr::pull(total_priors))

# ------------------------------------------------------------------------------
# Supplementary Figure 2
# ------------------------------------------------------------------------------
priors_vs_total_edges <- ggplot(tab, aes(x=total_edges, y=total_priors, color=group)) + 
  geom_hline(yintercept=0, color="grey") + 
  geom_point(alpha=0.5) + 
#  facet_wrap(~group, ncol=2, scales = "free_y") + 
  scale_color_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL)) + 
  scale_x_log10() + 
  background_grid(major="xy") +
  geom_smooth(method = "lm") + 
  labs(x="log10(Number of possible edges)",
       y="Total with prior information") +
   theme(legend.position = "none")
priors_vs_total_edges

priors_vs_total_edges_frac <- ggplot(tab, aes(x=total_edges, y=fraction_priors, color=group)) + 
  geom_hline(yintercept=0, color="grey") + 
  geom_point(alpha=0.5) + 
  #  facet_wrap(~group, ncol=2, scales = "free_y") + 
  scale_color_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL)) + 
  scale_x_log10() + 
  background_grid(major="xy") +
  geom_smooth(method = "lm") + 
  labs(x="log10(Number of possible edges)",
       y="Fraction with prior information") +
   theme(legend.position = "none")
priors_vs_total_edges_frac

prior_fraction <- ggplot(tab, aes(x=fraction_priors, fill=group)) +
  geom_histogram() + 
  scale_fill_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL)) + 
  facet_wrap(. ~ group, ncol=1) + 
  background_grid(major="xy") + 
  labs(x="Fraction of edges with priors")
prior_fraction

ga <- ggarrange(prior_fraction, 
                ggarrange(priors_vs_total_edges, priors_vs_total_edges_frac, 
                          ncol=1, labels=c("B", "C")),
                common.legend = T, legend="right",
                nrow=1, labels=c("A"))

save_plot("results/current/figures/prior_fractions.pdf", ga, 
          ncol=2, nrow=2)

prior_plot_thesis <- ggarrange(panel_d, 
                               ggarrange(priors_vs_total_edges,
                                         priors_vs_total_edges_frac,
                                         labels=c("B","C"), ncol=2), 
                               labels="A",
                               nrow = 2)

save_plot("results/current/figures/thesis/prior_plots.pdf", 
          prior_plot_thesis,
          ncol=2, nrow=2)

#-----------------------------------------------------------------------------
print("Prior plot using distinct prior categories (potential supp figure)")
# ------------------------------------------------------------------------------
get_prior_values <- function(files, type = c("snp-gene", "gene-gene", 
                                              "cpg-cpggene", "", "")) {
  
  unlist(sapply(files, function(prior_file) {
    
    # get priors and according ranges object
    pr <- readRDS(prior_file)
    ra <- readRDS(gsub("priors/", "ranges/", prior_file))
    
    # only meQTL have cpg-cpggene priors
    if(!grepl("meqtl", prior_file) & type == "cpg-cpggene") {
      return(NULL)
    }
    
    # get individual entities
    sentinel <- names(ra$sentinel)
    cpgs <- names(ra$cpgs)
    snp_genes <- ra$snp_genes$SYMBOL
    cpg_genes <- ra$cpg_genes$SYMBOL
    all_genes <- unique(c(ra$snp_genes$SYMBOL,
                   ra$cpg_genes$SYMBOL,
                   ra$trans_genes$SYMBOL))
    
    if(type == "snp-gene") {
      sub <- unlist(pr[sentinel, snp_genes])
    } else if(type == "gene-gene") {
      sub <- unlist(pr[all_genes[1], all_genes[2:length(all_genes)]])
    } else if(type == "cpg-cpggene") {
      sub <- unlist(pr[cpgs, cpg_genes])
    }
    
    # remove the pseudo priors
    out <- sub[sub > min(pr)]
    out
  }))
}

snp_gene <- tibble(prior=get_prior_values(finput, "snp-gene"), type="SNP-Gene")
gene_gene <- tibble(prior=get_prior_values(finput, "gene-gene"), type="Gene-Gene")
cpg_gene <- tibble(prior=get_prior_values(finput, "cpg-cpggene"), type="CpG-Gene")

data <- bind_rows(snp_gene, gene_gene, cpg_gene)

# plot densities
prior_densities <- data %>% ggplot(aes(x=prior, stat(density), col=type)) +
  geom_freqpoly(bins = 40) + 
  scb_priors + 
  scale_x_continuous(limits = c(0,1))

# ------------------------------------------------------------------------------
print("Figure 1 - Combine panels and saving.")
# ------------------------------------------------------------------------------
final_figure1 <- ggarrange(ggarrange(nullGrob(), panel_b1, panel_b2, ncol=3, align="v", 
                                  widths=c(0.02,0.49,0.49)),
                        ggarrange(panel_c, panel_d, ncol=2, align="v", 
                                  labels = c("C", "D")),
                        nrow=2, labels = c("B"))

save_plot("results/current/figures/figure1.pdf",
          plot=final_figure1, nrow = 2, ncol = 2, base_aspect_ratio = 1)

# ------------------------------------------------------------------------------
print("Figure 2 - Panel A")
# ------------------------------------------------------------------------------
# input directory containing individual simulation validation
# for the meQTL analysis
finput <- paste0(RESULT_PATH, "simulation/validation-subsetall.txt")

# create data-matrix
print("Reading simulation validation results...")
res <- read_tsv(finput) %>%
  mutate(R = paste0("R=", rdegree))


# create nicer method names
tab <- res %>% 
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("bdgraph_no_priors$", "bdgraph (empty)", comparison),
         comparison = gsub("bdgraph_no_priors_full$", "bdgraph (full)", comparison),
         comparison = gsub("glasso_no_priors","glasso", comparison)) %>%
  dplyr::rename(method=comparison)
toplot <- tab %>%
  filter(method != "bdgraph (full)") %>%
  mutate(is_prior_based = method %in% c("bdgraph (priors)", 
                                        "glasso (priors)",
                                        "irafnet")) %>%
  mutate(method = gsub("bdgraph \\(empty\\)", "bdgraph", method)) %>%
  mutate(method = factor(
    method,
    ordered = T,
    levels = c(
      "bdgraph",
      "bdgraph (priors)",
      "glasso",
      "glasso (priors)",
      "irafnet",
      "genie3",
      "genenet"
    )
  ))
  
# get the MCC plot
simulation_mcc <- ggplot(toplot,
                         aes(y=MCC, 
                             x=R, 
                             color=method)) +
  stat_boxplot(geom="errorbar", width=.75)+
  geom_boxplot(outlier.size=0, alpha=0.5, coef=0, outlier.shape = NA, ) + 
  stat_summary(fun.y=median, geom="smooth", 
               position=position_dodge(0.75),
               aes(group=method),lwd=0.8) +
  scb_graphs +
  #geom_boxplot(data = tab, alpha=0.5, aes(y=density_true, x=R), 
  #             fill="#666666", color="#666666",
  #             inherit.aes = FALSE,
  #             width=.3) +
  #scale_y_continuous(limits=c(min(tab$MCC),1), 
  #                   breaks=c(0,0.25,0.5,0.75,1),
  #                   sec.axis = sec_axis(trans = ~ ., 
  #                                       name="true graph density",
  #                                       breaks=seq(0,0.7,by=0.1))) +
  background_grid(major="x") +
  geom_vline(xintercept = 11.5, size=1.5, color="black", linetype="dashed") + 
  labs(x="prior error",
       y="MCC",
       fill="") + 
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14), 
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        axis.title.x = element_text(size=14, 
                                    margin = margin(-1, 0, 0, 0, unit = "lines")),
        plot.margin = margin(0.2,1,0.2,0.2,"cm"))
simulation_mcc

save_plot("results/current/figures/thesis/simulation_performance_thesis.pdf",
          simulation_mcc + theme(legend.position = "bottom"), base_asp = 2)

# ------------------------------------------------------------------------------
print("Method comparison based on reviewer comments (fixed error, subset size)")
toplot_filtered <- filter(toplot, R=="R=0.1")

# compare prior vs non_prior for glasso and bdgraph
toplot_glasso_bdgraph <- filter(toplot_filtered,
                                grepl("glasso|bdgraph", method)) %>%
  mutate(method_group = case_when(grepl("glasso", method) ~ "glasso",
                                  TRUE ~ "bdgraph"))

# compare overall prior VS non prior (should we do this?)
p1 <- ggplot(toplot_filtered %>%
               mutate(is_prior_based = as.factor(case_when(!is_prior_based ~ "no priors",
                                                           is_prior_based ~ "prior based"))),
             aes(
               y = MCC,
               x = is_prior_based,
               color = method,
               group = reorder(method, (-MCC))
             )) +
  #  geom_violin(draw_quantiles = c(0.5)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method.args = list(alternative = "less"),
                             comparisons = list(base = c(1, 2))) +
  scb_graphs + 
  theme(axis.title.x = element_blank())

p1

# ------------------------------------------------------------------------------
# check now the same as above, but progress over prior noise
toplot_filtered <- filter(toplot, R!="R=rbinom")

# compare prior vs non_prior for glasso and bdgraph
toplot_glasso_bdgraph <- filter(toplot_filtered,
                                grepl("glasso|bdgraph", method)) %>%
  mutate(method_group = case_when(grepl("glasso", method) ~ "glasso",
                                  TRUE ~ "bdgraph"))
p2 <- ggplot(toplot_glasso_bdgraph,
             aes(
               y = MCC,
               x = method,
               color = method,
               group = method
             )) +
  #  geom_violin(draw_quantiles = c(0.5)) +
  geom_boxplot() +
  facet_wrap(R~., strip.position = "left") +
  ggpubr::stat_compare_means(
    method.args = list(alternative = "two"),
    comparisons = list(bdgraph = c(1, 2),
                       glasso = c(3, 4))
  ) +
  scb_graphs +
  theme(axis.text.x = element_blank())

ga <- ggarrange(p2, ggarrange(p1, nullGrob(), ncol=2),  heights=c(0.7,0.3),
                nrow=2, labels = "AUTO", common.legend = F)
ga

save_plot("results/current/revisions/figures/prior_vs_non_prior_stat_sign_all_noise.pdf",
          ga, ncol=2, nrow = 3.3)
save_plot("results/current/revisions/figures/prior_vs_non_prior_stat_sign_all_noise.png",
          ga, ncol=2, nrow = 3.3)

# ------------------------------------------------------------------------------
print("Sample size simulation plot")

finput <- file.path(RESULT_PATH, "simulation/validation-subsets.txt")

tab <- read_tsv(finput) %>% 
  filter(rdegree == 0) %>%
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("bdgraph_no_priors$", "bdgraph (empty)", comparison),
         comparison = gsub("bdgraph_no_priors_full$", "bdgraph (full)", comparison),
         comparison = gsub("glasso_no_priors","glasso", comparison)) %>%
  dplyr::rename(method=comparison) %>%
  mutate(subset = factor(subset, level=c(seq(50,600,by=50)), ordered = T))

toplot <- tab %>%
  filter(method != "bdgraph (full)") %>%
  mutate(method = gsub("bdgraph \\(empty\\)", "bdgraph", method)) %>%
  mutate(method = factor(
    method,
    ordered = T,
    levels = c(
      "bdgraph",
      "bdgraph (priors)",
      "glasso",
      "glasso (priors)",
      "irafnet",
      "genie3",
      "genenet"
    )
  ))

simulation_sample_size_mcc <- ggplot(toplot,
                      aes(y = MCC,
                          x = subset,
                          color = method)) +
  scb_graphs +
  stat_boxplot(geom = "errorbar", width = .75) +
  geom_boxplot(
    outlier.size = 0,
    alpha = 0.5,
    coef = 0,
    outlier.shape = NA
  ) +
  stat_summary(
    fun.y = median,
    geom = "smooth",
    position = position_dodge(0.75),
    aes(group = method),
    lwd = 0.8
  ) +
  scale_y_continuous(
    limits = c(min(tab$MCC), 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
    
  ) +
  labs(color = "method",
       x = "subset size") + 
  theme(axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12))

simulation_sample_size_mcc
save_plot("results/current/figures/thesis/simulation_performance_sample_size_thesis.pdf",
          simulation_sample_size_mcc + theme(legend.position = "bottom"), base_asp = 2)

# ------------------------------------------------------------------------------
print("Figure 2 - Panel B")
# ------------------------------------------------------------------------------
# read the validation results for meqtls and eqtls and combine them
# TODO: change to rerun
meqtl_expr <- read_tsv(paste0(RESULT_PATH, "validation_expr/_rerun/validation_all_meqtl.txt"))
meqtl_tfa <- read_tsv(paste0(RESULT_PATH, "validation_tfa/_rerun/validation_all_meqtl.txt"))
meqtl <- bind_rows(meqtl_expr,meqtl_tfa) %>%
  mutate(type=c(rep("Expression", nrow(meqtl_expr)), rep("TF activities", nrow(meqtl_tfa))),
         qtl_type="meQTL")

eqtl_expr <- read_tsv(paste0(RESULT_PATH, "validation_expr/_rerun/validation_all_eqtlgen.txt"))
eqtl_tfa <- read_tsv(paste0(RESULT_PATH, "validation_tfa/_rerun/validation_all_eqtlgen.txt"))
eqtl <- bind_rows(eqtl_expr,eqtl_tfa) %>%
  mutate(type=c(rep("Expression", nrow(eqtl_expr)), rep("TF activities", nrow(eqtl_tfa))),
         qtl_type="eQTL")

data <- bind_rows(meqtl, eqtl)

data <- data %>% 
  mutate(graph_type = gsub("bdgraph$", "bdgraph (priors)", graph_type),
         graph_type = gsub("glasso$", "glasso (priors)", graph_type),
         graph_type = gsub("_no_priors", "", graph_type)) %>%
  dplyr::select(graph_type, type, graph_score, cross_cohort_mcc)

expr_plot <- data %>%
  filter(type %in% "Expression") %>%
  ggplot(aes(x=reorder(graph_type, -cross_cohort_mcc, FUN=median),
             y=cross_cohort_mcc)) +
  geom_boxplot() + 
  geom_point(position=position_jitter(width = 0.15),
             alpha=0.1) +
  background_grid(major="xy") +
  labs(title="",
       y="MCC",
       x="method") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        legend.position = "bottom",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))

#aes(ymin=cross_cohort_mcc-sd, ymax=cross_cohor_mcc+sd))

tfa_expr_plot <- data %>%
  ggplot(aes(x=reorder(graph_type, -cross_cohort_mcc, FUN=median), 
             y=cross_cohort_mcc, color=type)) + 
  #geom_violin(position = "dodge", draw_quantiles = 0.5, scale = "width") + 
  geom_boxplot(position="dodge") +
  #geom_point(position=position_jitterdodge(jitter.width = 0.15,
  #                                         dodge.width = 0.75),
  #           alpha=0.2) +
  scale_color_grey(start = 0.2, end = 0.7) + 
  #stat_summary(fun.y="median", geom = "line", aes(group=type), size=1)+
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  labs(title="",
       y="MCC",
       x="method",
       color = "measure:") + 
  background_grid(major="xy") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        legend.position = "right",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15), 
        plot.margin = margin(0.2,1,0.2,0.2,"cm")) + 
  stat_compare_means(aes(group=type), size = 10, method="wilcox.test", 
                     hide.ns=T, label = "p.signif", show.legend = F)
tfa_expr_plot
# possible alternatives?
data %>% ggplot(aes(x=cross_cohort_mcc)) + 
  geom_freqpoly(aes(color=type), stat = "density") +
  scb_binary +
  labs(color="measure:", x= "MCC") + background_grid(major="xy")
data %>% ggplot(aes(x=reorder(graph_type, -cross_cohort_mcc, median), 
                    y=cross_cohort_mcc, color=type)) + 
  geom_boxplot(position="dodge") + 
  stat_summary(fun.y="mean", geom = "line", aes(group=type))

save_plot("results/current/figures/thesis/tfa_vs_expr_thesis.pdf", 
          tfa_expr_plot + 
            scale_color_grey(start = 0.3, end = 0.6) + 
            theme(plot.margin = margin(0.6,1,0.2,0.2,"cm")), 
          base_asp = 2)

# This plot can be used as supplement
tfa_expr_graph_score <- data %>%
  ggplot(aes(x=reorder(graph_type, -graph_score, FUN=median), 
             y=graph_score, color=type)) + 
  #geom_violin(position = "dodge", draw_quantiles = 0.5, scale = "width") + 
  geom_boxplot(position="dodge", outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15,
                                           dodge.width = 0.75),
             alpha=0.2) +
  scale_color_grey(start = 0.3, end = 0.6) + 
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  labs(title="",
       y="graph score",
       x="method",
       color = "measure:") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        legend.position = c(0.7,0.8),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) + 
  stat_compare_means(aes(group=type), size=2, method="wilcox.test", 
                     hide.ns=T, label = "p.signif")

tfa_expr_graph_score
save_plot("results/current/supplement/graph_scores.pdf", tfa_expr_graph_score)


# ------------------------------------------------------------------------------
print("Figure 2 - Compile full plot")
# ------------------------------------------------------------------------------
figure2 <- ggarrange(simulation_mcc,
                     tfa_expr_plot, 
                     ncol = 2, labels = c("A", "B", "C"),
                     align="h")
figure2_alt <- ggarrange(ggarrange(simulation_sample_size_mcc, 
                                   simulation_mcc,  ncol=2, common.legend = T,
                                   legend="right", labels=c("A", "B"), align = "h"),
                         ggarrange(expr_plot,
                                   tfa_expr_plot,
                                   ncol=2, 
                                   common.legend = T, 
                                   legend="right", 
                                   labels=c("C","D")), 
                         nrow=2)
figure2
figure2_alt
save_plot("results/current/figures/figure2.pdf",
          plot=figure2, nrow = 2, ncol = 2, 
          base_aspect_ratio = 2.5)

save_plot("results/current/revisions/figures/figure2_alternative.pdf",
          plot=figure2_alt, nrow = 2, ncol = 2, 
          base_aspect_ratio = 2)

# ------------------------------------------------------------------------------
print("SF1 - graph scores")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("SFX - number of inferred edges")
# ------------------------------------------------------------------------------
meqtl <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_meqtl.txt") %>% mutate(type = "meqtl")
eqtl <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_eqtlgen.txt") %>% mutate(type = "eqtl")
result <- bind_rows(eqtl, meqtl)
result %>% ggplot(aes(x=number_edges, color=graph_type)) + 
  geom_freqpoly() + 
  facet_wrap(.~cohort) + 
  scale_x_log10() + 
  background_grid(major="xy")
result %>% 
  filter(graph_type %in% c("irafnet", "genie3")) %>% 
  ggplot(aes(x=number_edges, color=graph_type)) + 
  geom_freqpoly() + 
  facet_wrap(.~cohort) + 
  background_grid(major="xy") + 
  scale_x_log10()

# ------------------------------------------------------------------------------
print("SFXX - SNp-gene simulation evaluation")
# ------------------------------------------------------------------------------
# read in summary tables
finput <- list.files("results/current/biogrid_stringent/simulation_rerun/validation/snp_gene_recovery", 
                     "*.tsv",
                     full.names = T)
res <- lapply(finput, function(f) { read_tsv(f)}) %>% 
  bind_rows() %>%
  mutate(graph_type = gsub("bdgraph$", "bdgraph (priors)", graph_type),
         graph_type = gsub("glasso$", "glasso (priors)", graph_type),
         graph_type = gsub("bdgraph_no_priors$", "bdgraph (empty)", graph_type),
         graph_type = gsub("bdgraph_no_priors_full$", "bdgraph (full)", graph_type),
         graph_type = gsub("glasso_no_priors","glasso", graph_type)) %>%
  filter(graph_type != "bdgraph (full)")

gp1 <- res %>% 
  ggplot(aes(y = ACC, x= reorder(graph_type, -ACC, median))) + 
  geom_boxplot() + 
  labs(x = "") + 
  theme(axis.text.x = element_text(angle = -45, hjust=0))

# get summary for each graph_type (total number of snp genes in true graphs and
# identified)
gp2 <- res %>%
  group_by(graph_type) %>%
  summarise(total_snp_genes = sum(number_snp_genes),
            total_recovered = sum(recovered_snp_genes)) %>%
  mutate(total_fraction = total_recovered / total_snp_genes) %>%
  ggplot(aes(x=reorder(graph_type, -total_fraction, max), y = total_fraction)) + 
  geom_bar(stat="identity") +
  scale_y_continuous(limits=c(0,1)) + 
  labs(x = "", y = "overall fraction recovered") + 
  theme(axis.text.x = element_text(angle = -45, hjust=0), 
        plot.margin = margin(0.5,1,0.5,0.5,"cm"))

save_plot("results/current/figures/thesis/snp_gene_recovery_thesis.pdf", gp2)

final <- ggpubr::ggarrange(gp1,gp2, align = "h", ncol = 2, labels = "AUTO")
save_plot("results/current/supplement/snp_gene_recovery.pdf", final, 
          ncol=2, 
          base_asp = 1.2)
# ------------------------------------------------------------------------------
print("Done.\nSessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()





# ------------------------------------------------------------------------------
print("Matthias evaluation")
# ------------------------------------------------------------------------------
finput <- paste0(RESULT_PATH, "simulation_rerun/validation-subsetall.txt")

# create data-matrix
print("Reading simulation validation results...")
res <- read_tsv(finput) %>%
  mutate(R = paste0("R=", rdegree))


# create nicer method names
tab <- res %>% 
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("bdgraph_no_priors$", "bdgraph (empty)", comparison),
         comparison = gsub("bdgraph_no_priors_full$", "bdgraph (full)", comparison),
         comparison = gsub("glasso_no_priors","glasso", comparison)) %>%
  dplyr::rename(method=comparison)
toplot <- tab %>%
  filter(method != "bdgraph (full)") %>%
  mutate(method = gsub("bdgraph \\(empty\\)", "bdgraph", method)) %>%
  mutate(method = factor(
    method,
    ordered = T,
    levels = c(
      "bdgraph",
      "bdgraph (priors)",
      "glasso",
      "glasso (priors)",
      "irafnet",
      "genie3",
      "genenet"
    )
  ))

toplot_all_no_priors <- toplot %>%
  filter(method %in% c("bdgraph", "glasso", "genie3", "genenet"))

toplot_priors <- toplot %>%
  filter(method %in% c("bdgraph","glasso","bdgraph (priors)", "glasso (priors)"))

simulation_mcc_no_priors <- ggplot(toplot_all_no_priors,
                         aes(y=MCC, 
                             x=R, 
                             color=method)) +
  stat_boxplot(geom="errorbar", width=.75)+
  geom_boxplot(outlier.size=0, alpha=0.5, coef=0, outlier.shape = NA, ) + 
  stat_summary(fun.y=median, geom="smooth", 
               position=position_dodge(0.75),
               aes(group=method),lwd=0.8) +
  scb_graphs +
  background_grid(major="x") +
  geom_vline(xintercept = 11.5, size=1.5, color="black", linetype="dashed") + 
  labs(x="prior error",
       y="MCC",
       fill="") + 
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14), 
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        axis.title.x = element_text(size=14, 
                                    margin = margin(-1, 0, 0, 0, unit = "lines")),
        plot.margin = margin(0.2,1,0.2,0.2,"cm"))

simulation_mcc_priors <- ggplot(toplot_priors,
                                   aes(y=MCC, 
                                       x=R, 
                                       color=method)) +
  stat_boxplot(geom="errorbar", width=.75)+
  geom_boxplot(outlier.size=0, alpha=0.5, coef=0, outlier.shape = NA, ) + 
  stat_summary(fun.y=median, geom="smooth", 
               position=position_dodge(0.75),
               aes(group=method),lwd=0.8) +
  scb_graphs +
  background_grid(major="x") +
  geom_vline(xintercept = 11.5, size=1.5, color="black", linetype="dashed") + 
  labs(x="prior error",
       y="MCC",
       fill="") + 
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14), 
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12),
        axis.title.x = element_text(size=14, 
                                    margin = margin(-1, 0, 0, 0, unit = "lines")),
        plot.margin = margin(0.2,1,0.2,0.2,"cm"))

save_plot("results/current/figures/simulation_performance_matthias_eval_no_priors.pdf",
          simulation_mcc_no_priors + theme(legend.position = "bottom"), base_asp = 2)
save_plot("results/current/figures/simulation_performance_matthias_eval_priors.pdf",
          simulation_mcc_priors + theme(legend.position = "bottom"), base_asp = 2)
# ------------------------------------------------------------------------------
print("Sample size simulation plot")

finput <- file.path(RESULT_PATH, "simulation_rerun/validation-subsets.tsv")

tab <- read_tsv(finput) %>% 
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("bdgraph_no_priors$", "bdgraph (empty)", comparison),
         comparison = gsub("bdgraph_no_priors_full$", "bdgraph (full)", comparison),
         comparison = gsub("glasso_no_priors","glasso", comparison)) %>%
  dplyr::rename(method=comparison) %>%
  mutate(subset = factor(subset, level=c(seq(50,600,by=50)), ordered = T))

toplot <- tab %>%
  filter(method != "bdgraph (full)") %>%
  mutate(method = gsub("bdgraph \\(empty\\)", "bdgraph", method)) %>%
  mutate(method = factor(
    method,
    ordered = T,
    levels = c(
      "bdgraph",
      "bdgraph (priors)",
      "glasso",
      "glasso (priors)",
      "irafnet",
      "genie3",
      "genenet"
    )
  ))

toplot_all_no_priors <- toplot %>%
  filter(method %in% c("bdgraph", "glasso", "genie3", "genenet"))

toplot_priors <- toplot %>%
  filter(method %in% c("bdgraph","glasso","bdgraph (priors)", "glasso (priors)"))

simulation_sample_size_mcc_no_priors <- ggplot(toplot_all_no_priors %>% filter(rdegree==0),
                                     aes(y = MCC,
                                         x = subset,
                                         color = method)) +
  scb_graphs +
  stat_boxplot(geom = "errorbar", width = .75) +
  geom_boxplot(
    outlier.size = 0,
    alpha = 0.5,
    coef = 0,
    outlier.shape = NA
  ) +
  stat_summary(
    fun.y = median,
    geom = "smooth",
    position = position_dodge(0.75),
    aes(group = method),
    lwd = 0.8
  ) +
  scale_y_continuous(
    limits = c(min(tab$MCC), 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
    
  ) +
  labs(color = "method",
       x = "subset size") + 
  theme(axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12))

simulation_sample_size_mcc_priors <- ggplot(toplot_priors %>% filter(rdegree==0),
                                               aes(y = MCC,
                                                   x = subset,
                                                   color = method)) +
  scb_graphs +
  stat_boxplot(geom = "errorbar", width = .75) +
  geom_boxplot(
    outlier.size = 0,
    alpha = 0.5,
    coef = 0,
    outlier.shape = NA
  ) +
  stat_summary(
    fun.y = median,
    geom = "smooth",
    position = position_dodge(0.75),
    aes(group = method),
    lwd = 0.8
  ) +
  scale_y_continuous(
    limits = c(min(tab$MCC), 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
    
  ) +
  labs(color = "method",
       x = "subset size") + 
  theme(axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45, size=12))


save_plot("results/current/figures/simulation_performance_samplesize_matthias_eval_no_priors.pdf",
          simulation_sample_size_mcc_no_priors + theme(legend.position = "bottom"), base_asp = 2)
save_plot("results/current/figures/simulation_performance_samplesize_matthias_eval_priors.pdf",
          simulation_sample_size_mcc_priors + theme(legend.position = "bottom"), base_asp = 2)
