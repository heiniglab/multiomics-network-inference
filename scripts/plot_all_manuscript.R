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
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
# get the hg19 chromosome definitions
library(BSgenome.Hsapiens.UCSC.hg19)
hg19info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

# set up theme and colors
theme_set(theme_cowplot())
theme_update(legend.text = element_text(size=10), 
             legend.title=element_text(size=10))

sfb_graphs <- scale_fill_brewer(palette="Set2")
scb_graphs <- scale_color_brewer(palette="Set2")
sfb_binary <- scale_fill_brewer(palette = "Accent")
scb_binary <- scale_color_brewer(palette = "Accent")
bgm <- background_grid(major = "xy")
group_cols <- brewer.pal("Set2", n=3)
COLORS <- list(MEQTL = group_cols[1],
               EQTL = group_cols[2])

# ------------------------------------------------------------------------------
print("Figure 1 - Panel A")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Now we create the plot showing the overview over the locus graphs
# ------------------------------------------------------------------------------
# get meqtl ranges
fmeqtl_ranges <- list.files("results/current/biogrid_stringent/ranges/", "*_meqtl.rds", 
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
feqtl_ranges <- list.files("results/current/biogrid_stringent/ranges/", "*_eqtlgen.rds", 
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

# ------------------------------------------------------------------------------
# Plot the distributions
# ------------------------------------------------------------------------------
toplot <- melt(counts) %>% 
  arrange(value)

panel_b <- ggplot(toplot, aes(color=group, x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=23, position = position_jitterdodge(jitter.width = 0.15)) + 
  scale_color_manual(values=c(meQTL = COLORS$MEQTL, eQTL = COLORS$EQTL, TWAS = COLORS$TWAS)) + 
  scale_y_log10() +
  xlab("entity type") + 
  ylab("count") + theme(plot.margin = margin(1,0.2,0.2,0.2, unit="lines"))

# ------------------------------------------------------------------------------
print("Figure 1 - Panel A")
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
    geom_point(aes(x=x_bin, y=count), color=color, shape=23) + 
    scale_x_continuous(expand = c(0.01, 0.01), breaks = as.vector(breaks), labels = NULL) + 
    xlab("") + 
    background_grid(major = "xy")
  
  ymp <- ggplot(y_margin) + 
    geom_point(aes(x=y_bin, y=count), color=color, shape=23) + 
    coord_flip() + scale_x_continuous(
      expand = c(0.01, 0.01),
      breaks = as.vector(breaks),
      labels = NULL) + 
    xlab("") + 
    theme(axis.text.x = element_text(angle=90, 
                                     vjust=0.5, 
                                     hjust=1)) + 
    background_grid(major = "xy")
  
  # the discretized plot
  gd <-
    ggplot(pairs_binned_discrete) + 
    geom_tile(..., fill=color,
              aes(x = x_bin, y = y_bin)) +
    theme(
      text = element_text(size = 11),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(vjust = 0.5, angle = 90),
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

tile_width <- 10

eqtl_plot <- plot_pairs(hg19info, 
                        eqtl_regions, 
                        eqtl_regions$trans_ranges, 
                        log = F, "SNPs", "trans Genes", resolution = 1e6, width=tile_width, color=COLORS$EQTL)
meqtl_plot <- plot_pairs(hg19info, 
                         meqtl_regions, 
                         meqtl_regions$trans_ranges, 
                         log=F, "SNPs", "CpGs", resolution = 1e6, width=tile_width, color=COLORS$MEQTL)

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
#meqtl_ggplot <- ggdraw() + 
#  draw_plot(arrangeGrob(nullGrob(), 
#                        plot_grid(bp+theme(plot.background = element_rect(fill="transparent")), 
#                                  cp, ncol=2, align="h", rel_widths=wr), 
#                        nrow=2, heights = hr)) +
#  draw_plot(arrangeGrob(plot_grid(ap,bp, nrow=2, align="v", rel_heights=hr), 
#                        nullGrob(), ncol=2, widths = wr))

# possible alternative:
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
#eqtl_ggplot <- ggdraw() + 
#  draw_plot(arrangeGrob(nullGrob(), 
#                        plot_grid(bp+theme(plot.background = element_rect(fill="transparent")), cp, ncol=2, align="h", rel_widths=wr), 
#                        nrow=2, heights = hr)) +
#  draw_plot(arrangeGrob(plot_grid(ap,bp, nrow=2, align="v", rel_heights=hr), 
#                        nullGrob(), ncol=2, widths = wr))

eqtl_ggplot <- ggarrange(ggarrange(ap, nullGrob(), widths = wr, ncol=2), 
                          ggarrange(bp, cp, widths=wr, ncol=2, align="h"), 
                          nrow=2, heights = hr, align="v")

panel_a <- grid.arrange(meqtl_ggplot, 
                        eqtl_ggplot,
                        ncol=2)

# extract legend for wide plot
panel_b <- panel_b + 
  theme(legend.text=element_text(size=rel(1)),
        legend.title=element_text(size=rel(1.5)),
        legend.key.size = unit(2,"lines"))

legend <- cowplot::get_legend(panel_b)
panel_b <- panel_b + theme(legend.position = "none")

final_wide <- ggarrange(ggarrange(panel_a, panel_b, nrow=2), legend, 
                        ncol=2, widths = c(0.8,0.2))
final_wide

save_plot("figure1_wide.pdf",
          plot=final_wide, nrow = 2, base_aspect_ratio = 4)

# ------------------------------------------------------------------------------
print("Figure 2 - Panel A")
# ------------------------------------------------------------------------------
# input directory containing individual simulation validation
# for the meQTL analysis
dinput <- "results/current/biogrid_stringent/simulation/validation/"
finput <- list.files(dinput, "*.txt", full.names = T)

# create data-matrix
temp <- lapply(finput, function(f) {
  f <- read_tsv(f)
  if(nrow(f) > 0) {
    f <- f %>%
      mutate(R = paste0("R=", rdegree)) %>%
      rename(name = X1)
  }
})
tab <- bind_rows(temp)

# create nicer method names
tab <- tab %>% 
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("_no_priors", "", comparison))

# get the MCC plot
simulation_mcc <- ggplot(tab,
                         aes(y=MCC, 
                             x=R, 
                             color=reorder(comparison, -MCC, median))) +
  stat_boxplot(geom="errorbar", width=.75)+
  geom_boxplot(outlier.size=0, alpha=0.5, coef=0, outlier.shape = NA) + 
  stat_summary(fun.y=median, geom="smooth", 
               position=position_dodge(0.75),
               aes(group=comparison),lwd=0.8) +
  scb_graphs +
  geom_boxplot(data = tab, alpha=0.5, aes(y=density_true, x=R), 
               fill="#666666", color="#666666",
               inherit.aes = FALSE,
               width=.3) +
  scale_y_continuous(limits=c(min(tab$MCC),1), 
                     sec.axis = sec_axis(trans = ~ ., 
                                         name="true graph density",
                                         breaks=seq(0,0.7,by=0.1))) +
  background_grid(major="xy") +
  labs(x="prior error",
       y="MCC",
       fill="", color="method") + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
print("Figure 2 - Panel B")
# ------------------------------------------------------------------------------
# read the validation results for meqtls and eqtls and combine them
meqtl_expr <- read_tsv("results/current/biogrid_stringent/validation_expr/validation_all_meqtl.txt")
meqtl_tfa <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_meqtl.txt")
meqtl <- bind_rows(meqtl_expr,meqtl_tfa) %>%
  mutate(type=c(rep("expr", nrow(meqtl_expr)), rep("tfa", nrow(meqtl_tfa))),
         qtl_type="meQTL")

eqtl_expr <- read_tsv("results/current/biogrid_stringent/validation_expr/validation_all_eqtlgen.txt")
eqtl_tfa <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_eqtlgen.txt")
eqtl <- bind_rows(eqtl_expr,eqtl_tfa) %>%
  mutate(type=c(rep("expr", nrow(eqtl_expr)), rep("tfa", nrow(eqtl_tfa))),
         qtl_type="eQTL")

data <- bind_rows(meqtl, eqtl)

tfa_expr_plot <- data %>%
  ggplot(aes(x=reorder(graph_type, -cross_cohort_mcc, FUN=median), 
             y=cross_cohort_mcc, color=type)) + 
  geom_boxplot(position = "dodge") + 
  scb_binary + 
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  labs(title="",
       y="MCC",
       x="method",
       color = "measure:") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45),
        legend.position = "bottom")
tfa_expr_plot

# ------------------------------------------------------------------------------
print("Figure 2 - Compile full plot")
# ------------------------------------------------------------------------------
ggarrange(simulation_mcc,
          tfa_expr_plot,
          ncol = 2, labels = c("A", "B"))

# ------------------------------------------------------------------------------
print("Done.\nSessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
