#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
#' 
#' Script gathering all validation results and plotting them
#'
#' @author Johann Hawe johann.hawe@helmholtz-muenchen.de
#'
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load libraries and source scripts
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(scales)
library(knitr)
library(gridExtra)
library(graph)
source("scripts/lib.R")
cols <- set_defaultcolors()

# prepare some ggplot stuff
sfm <- scale_fill_manual(values=cols)
theme_set(theme_bw())

# ------------------------------------------------------------------------------
# Get snakemake params if available 
# ------------------------------------------------------------------------------

# inputs
fsummary <- snakemake@input[[1]]

# outputs
fstats <- snakemake@output$stats
fcratios <- snakemake@output$cratios
fexpr <- snakemake@output$expr
fgene_types <- snakemake@output$gene_types
fmediation <- snakemake@output$mediation
fmediation_perc <- snakemake@output$mediation_perc
fmediation_distr <- snakemake@output$mediation_distr
fperf <- snakemake@output$perf

# ------------------------------------------------------------------------------
# load the large result table for all loci
# TODO we currently remove the results for the GO enrichment
# ------------------------------------------------------------------------------
tab <- read.table(fsummary, header=T, sep="\t", stringsAsFactors=F)

# ------------------------------------------------------------------------------
# create some basic plots of the gene counts
# ------------------------------------------------------------------------------
toplot <- tab[,c("sentinel", "cohort", "graph_type",
                 "number_nodes", "number_edges", 
                 "graph_density", "cluster", "cluster_sizes", 
                 "snp_cluster",
                 "snp_genes", "snp_genes_selected", 
                 "cpg_genes", "cpg_genes_selected", 
                 "tfs", "tfs_selected",
                 "spath", "spath_selected")]
toplot$snp_in_network <- !is.na(toplot[,"snp_cluster"])

# whether the SNP has been selected or not
gp1 <- ggplotGrob(ggplot(data=toplot, aes(snp_in_network)) + geom_histogram(stat="count") + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Number of networks in which the SNP has been selected at all"))

# show the distribution of number of nodes per network
gp2 <- ggplotGrob(ggplot(data=toplot, aes(number_nodes)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Distribution of the number of nodes in the networks"))

# show the distribution of number of edges per network
gp3 <- ggplotGrob(ggplot(data=toplot, aes(number_edges)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Distribution of the number of edges in the graph"))

# show the distribution of graph_densities network
gp4 <- ggplotGrob(ggplot(data=toplot, aes(graph_density)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Distribution of graph densities over all networks", "density= 2*|E| / |V|*(|V|-1)"))

# show the distribution of resulting clusters (number of clusters, largest cluster size)
cluster_ratios <- c()
for(i in 1:nrow(toplot)) {
  # get cluster ratio
  cluster_sizes <- toplot[i,"cluster_sizes"]
  largest_cluster <- sort(as.numeric(strsplit(cluster_sizes, ",")[[1]]), decreasing = T)[1]
  cluster_ratios <- c(cluster_ratios, 
                      largest_cluster/toplot[i,"number_nodes"])
}
toplot$cluster_ratio <- cluster_ratios

# summarize the first four plots in one file
ggsave(plot=grid.arrange(gp1, gp2, gp3, gp4, ncol=2),
       file=fstats, 
       width=12, 
       height=8)

# save the cluster ratios plot individually
gp <- ggplotGrob(ggplot(data=toplot, aes(cluster_ratio)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Ratio of the amount of nodes in the largest clusters vs all nodes in the network"))
ggsave(plot=gp,
       file=fcratios,
       width=10, height=10)

# ------------------------------------------------------------------------------
# Check how many genes are retained in the graphs for each gene_type
# ------------------------------------------------------------------------------

# get ratios
toplot$snp_gene_ratio <- toplot$snp_genes_selected / toplot$snp_genes
toplot$cpg_gene_ratio <- toplot$cpg_genes_selected / toplot$cpg_genes
toplot$tf_ratio <- toplot$tfs_selected / toplot$tfs
toplot$spath_ratio <- toplot$spath_selected / toplot$spath

# ------------------------------------------------------------------------------
# Plot everything
# ------------------------------------------------------------------------------

# snp gene ratio
use <- toplot$snp_genes_selected>0
ggp1 <- ggplotGrob(ggplot(data=toplot[use,,drop=F], aes(snp_gene_ratio)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Ratio of number of selected SNP-genes vs all SNP-genes"))
# cpg gene ratio
use <- toplot$cpg_genes_selected>0
ggp2 <- ggplotGrob(ggplot(data=toplot[use,,drop=F], aes(cpg_gene_ratio)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Ratio of number of selected CpG-genes vs all CpG-genes"))
# tf ratio
use <- toplot$tfs_selected>0
ggp3 <- ggplotGrob(ggplot(data=toplot[use,,drop=F], aes(tf_ratio)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) + sfm +
  ggtitle("Ratio of number of selected TFs vs all TFs"))
# spath ratio
use <- toplot$spath_selected>0
ggp4 <- ggplotGrob(ggplot(data=toplot[use,,drop=F], aes(spath_ratio)) + geom_histogram() + 
  facet_grid(cohort ~ graph_type) +sfm +
  ggtitle("Ratio of number of selected shortest path genes vs all shortest path genes"))

# arrange and plot
ggsave(plot=grid.arrange(ggp1, ggp2, ggp3, ggp4, ncol=2),
       file=fgene_types,
       width=12, height=8)

#Above we see a simple overview over the selected entities (snp genes, cpg genes,
#transcription factors and shortest path genes).
#Shown is the log10-ratio of the number of entities selected by the GGM and the total 
#number of entities in the network. Therefore, low values indicate that comparatively more genes
#were dropped during the inference process and only fewer were selected for the final network.

#Below the mediation results on the different cohorts for all hotspots 
#in which at least one SNP gene was extracted by our algorithm are shown.
#The *mediation_selected* group shows the significance of the mediation analysis of 
#those selected genes, the *mediation* group shows the significance of all NOT selected genes. 
#Note, that for the 'mediation_selected' group, we always take the largest obtained p-value, whereas
#for the NOT selected group we always select the lowest obtained p-value as of now.

# ------------------------------------------------------------------------------
# Mediation summary plots
# ------------------------------------------------------------------------------
mediation <- tab[,c("sentinel", "cohort", "graph_type",
                    "mediation_notselected", "mediation_selected","log10_mediation")]
# for now we ignore results where we didn't have SNP genes at all
mediation <- mediation[!is.infinite(mediation$log10_mediation), ,drop=F]

# plot the mediation results
toplot <- mediation[order(mediation$mediation_selected),]
toplot$sentinel <- factor(toplot$sentinel, levels=unique(toplot$sentinel))
toplot <- melt(toplot, measure.vars=c(4,5), value.name = "pval")

gp1 <- ggplot(data=toplot, aes(x=graph_type, y=-log10(pval), fill=variable)) +
  geom_boxplot(outlier.color = NA) + 
  facet_grid(cohort ~ .) + sfm +
  geom_hline(yintercept = -log10(0.05), color="red") + 
  ggtitle("Mediation results of all SNP genes over all loci.")

# in addition we plot the log fold-change
# sort by log_fc
toplot <- mediation[order(mediation$log10_mediation, decreasing = T),]
toplot$sentinel <- factor(toplot$sentinel, levels=unique(toplot$sentinel))

toplot <- melt(toplot, measure.vars=6, value.name="log_fc")
toplot$favored <- factor(sign(toplot$log_fc),labels = c("not selected", "selected"))
gp2 <- ggplot(toplot, aes(x=sentinel, y=log_fc, fill=favored)) +
  geom_bar(stat="identity") + sfm +
  facet_grid(cohort ~ graph_type) +
  theme(axis.text.x = 
          element_text(angle = 90, hjust = 1, vjust=0)) +
  ylab("log10(all / selected)") + 
  ggtitle("Mediation results of all SNP genes, log-foldchanges.")

# save the two mediation plots
gp1 <- ggplotGrob(gp1)
gp2 <- ggplotGrob(gp2)

ggsave(plot=grid.arrange(gp1, gp2, nrow=2),
       file=fmediation,
       width=8, height=12)

#This figure shows a different summary of the mediation results. Shown is the log10 fold change
#of the mediation p-value for the not selected SNP genes ($p_n$) over the p-value for selected
#genes ($p_s$). The red bars (*not selected*) indicate fold changes where $p_n$ is lower 
#than the corresponding $p_s$, whereas the blue bars indicate the opposite. Overall, 
#**`r format(perc*100,digits=4)`\%** of all fold changes show negative fold changes 
#(i.e. $p_s$ being smaller than their respective $p_n$s). Currently, we select the
#minimal p-value obtained from the genes not selected via a GGM and the maximal p-value
#from the selected genes which is a rather conservative estimate. We could think about
#changing the way of summarizing the p-values...

# ------------------------------------------------------------------------------
# Percentage of mediating genes
# Here we look at the percentage of significant mediation results of the GGM selected
# SNP genes versus total number of SNP genes.
# ------------------------------------------------------------------------------

perc <- table(toplot$favored, toplot$graph_type)
perc <- perc["selected",] / colSums(perc)
print("Percentage of selected vs not selected mediating genes:")
print(perc)

# create data frame with all needed values
results <- tab[,c("graph_type", "cohort", "snp_genes", "snp_genes_selected",
		  "mediation_total",
                  "mediation_notselected_significant", "mediation_selected_significant")]

# get SNP gene percentages regarding mediation
results$sign_selected <- results[,"mediation_selected_significant"] / results[,"snp_genes_selected"]
results$snp_genes_notselected <- results[,"snp_genes"] - results[,"snp_genes_selected"]
results$sign_notselected <- results[,"mediation_notselected_significant"] / results[,"snp_genes_notselected"]

fg <- facet_grid(. ~ cohort)
sc <- scale_y_continuous(limits=c(0,1))
gg <- ggplot(data=results, aes(y=sign_selected, x=graph_type, fill=graph_type))
# plot distribution  of selected SNP genes showing mediation
gp1 <- gg + geom_violin(draw_quantiles=c(.25,.5,.75)) + 
			sfm + fg + sc +
			ggtitle("Percentage of selected SNP-genes showing mediation.")

# plot distribution  of not-selected SNP genes showing mediation
gp2 <- gg + aes(y=sign_notselected) +
       geom_violin(draw_quantiles=c(.25,.5,.75)) +
       sfm + fg + sc +
       ggtitle("Percentage of not-selected SNP-genes showing mediation.")

# plot distribution of percentage of total number of mediating genes
gp3 <- gg + aes(y=mediation_total/snp_genes) +
       geom_violin(draw_quantiles=c(.25,.5,.75)) +
       sfm + fg + sc +
       ggtitle("Percentage of SNP-genes showing mediation.")

# plot distribution  of total number of mediating genes
gp4 <- gg + aes(y=mediation_total) +
       geom_violin(draw_quantiles=c(.25,.5,.75)) +
       sfm + fg +
       ggtitle("Total number of SNP-genes showing mediation.")

ga <- grid.arrange(gp1,gp2,gp3,gp4,ncol=2)
ggsave(ga, file=fmediation_distr, width=12, height=9)

# handle different graph types separately
results <- split(results, f = results$graph_type)

temp <- lapply(names(results), function(n) {
  rsub <- results[[n]]
  cohort <- rsub[,2]
  
  summary <- lapply(unique(cohort), function(co) {
    # summarize results per cohort
    r <- colSums(rsub[rsub$cohort==co,c(-1,-2)], na.rm=T)
  
    # proportion of significant genes per group
    sign.selected <- r["mediation_selected_significant"] / r["snp_genes_selected"]
    sign.notselected <- r["mediation_notselected_significant"] / r["snp_genes_notselected"]
    
    # number of  
    toplot <- data.frame(selected=c(1,sign.selected), 
                         not.selected=c(1,sign.notselected),
			 cohort=c(co,co))
    rownames(toplot) <- c("total tested", "proportion significant")
  
    toplot <- melt(cbind(toplot, mediation = rownames(toplot)), 
                   id.vars = c('mediation', 'cohort'))
    return(toplot)
  })
  toplot <- do.call(rbind, summary)
  p <- ggplot(toplot, aes(x=variable,y=value, fill=mediation)) + 
    geom_bar(position="identity", stat="identity") +
    facet_grid(. ~ cohort) + sfm +
    theme(axis.text.x = 
             element_text(angle = 90, hjust = 1, vjust=0)) +
    ylab("percentage of genes") +
    xlab("mediation group") +
    scale_y_continuous(labels=percent_format()) + 
    ggtitle(n)
  return(p)
})
theme_update(plot.title = element_text(hjust = 0.5))
ga <- grid.arrange(temp[[1]], temp[[2]], temp[[3]], ncol=2)
ggsave(plot=ga, file=fmediation_perc,
       width=12, height=8)

# ------------------------------------------------------------------------------
# Specificity and sensitivity of gene selection over all sentinels
# Here we look at all sentinels and summarize the selected/not selected SNP 
# genes and their mediation values by calculating all TPs, TNs, FPs, and FNs.
# We do this for the individual cohorts as well as summed up over both
# ------------------------------------------------------------------------------
cohorts <- c("kora", "lolipop", "both")
graph_types <- unique(tab$graph_type)
values <- lapply(cohorts, function(cohort) {
  temp <- lapply(graph_types, function(graph_type) {
    if(cohort=="both") {
      cohort2 <- c("kora", "lolipop")
    } else { 
      cohort2 <- cohort
    }
    dat <- tab[tab$cohort %in% cohort2 & tab$graph_type == graph_type,
               c("sentinel", "snp_genes", "snp_genes_selected",
                 "mediation_notselected_significant", "mediation_selected_significant")]
    result <- sapply(1:nrow(dat), function(i) {
      d <- dat[i,,drop=F]
      v1 <- d$mediation_notselected_significant==0 & d$snp_genes_selected == 0
      v2 <- d$mediation_notselected_significant==0 & d$snp_genes_selected > 0
      v3 <- d$mediation_notselected_significant>0 & d$snp_genes_selected == 0 
      v4 <- d$mediation_selected_significant > 0
      # todo
      v5 <- ""
      return(c(v1,v2,v3,v4))  
    })
  
    result <- matrix(unlist(result), byrow = T, ncol=4)
    result <- colSums(result)
    names(result) <- c("TN", "FP", "FN", "TP")
    result
  })
  names(temp) <- paste(cohort, graph_types, sep=".")
  do.call(rbind, temp)
})
names(values) <- cohorts
values <- as.data.frame(do.call(rbind, values))
values$cohort <- unlist(lapply(strsplit(rownames(values), "\\."), "[[", 1))
values$graph_type <- unlist(lapply(strsplit(rownames(values), "\\."), "[[", 2))

# define performance summary methods
sens <- function(d) { d[,"TP"] / (d[,"TP"] + d[,"FN"]) }
spec <- function(d) { d[,"TN"] / (d[,"TN"] + d[,"FP"]) }
prec <- function(d) { d[,"TP"] / (d[,"TP"] + d[,"FP"]) }
f1 <- function(d) { 2 * 1 / ((1/sens(d)) + 
			     (1/prec(d))) }
# ------------------------------------------------------------------------------
# save a simple plot showing the performance values
# ------------------------------------------------------------------------------
df <- values
df$sensitivity <- sens(df)
df$specificity <- spec(df)
df$f1 <- f1(df)

df <- melt(df, measure.vars=c(7,8,9), value.name="performance")
gp <- ggplot(data=df, aes(y=performance, x=graph_type, fill=graph_type)) + 
	geom_bar(stat="identity") + 
	facet_grid(cohort ~ variable) + scale_y_continuous(limits=c(0,1)) +
        sfm +
	ggtitle("Performance of GGMs on different cohorts", "Baseline defined via significant mediation genes. Summary over all sentinels.")

ggsave(plot=gp, file=fperf,
       width=10,height=12)

# ------------------------------------------------------------------------------
# CpG-Gene validation
# ------------------------------------------------------------------------------
#TODO
#hist(tab$bonder_cis_eQTM)

# ------------------------------------------------------------------------------
# Gene-Gene validation 
#
#The gene-gene links were validated using external expression data from GEO (ARCHS4),
#Geuvadis as well as the data from either LOLIPOP or KORA, depending on which cohort 
#the bGGM was calculated.

#To assess how well our approach recovers co-expression of genes in the independent
#datasets, we calculate all gene-vs-gene (gvg) correlations in those data and create a
#confusion matrix as follows:

#temp <- matrix(c("W","X","Y","Z"),nrow=2,ncol=2)
#colnames(temp) <- c("In GGM", "Not in GGM")
#rownames(temp) <- c("Correlation sign.","Correlation not sign.")
#kable(temp,align = c("c","c"))

#We deem a correlation between two genes to be significant if 1) $p-value < 0.01$
#and 2) $\rho > 0.3$ ($\rho$ being the Pearson Correlation Coefficient).

#We test whether we can reject $H0$ (i.e. whether we detect a correlation regardless of it 
#being significant in the independent data) using a Fisher's exact test.
#Shown below are the p-values for those tests over all loci on the different datasets.
#The black lines indicate a significance level of $0.05$.
# ------------------------------------------------------------------------------

# get the relevant data (pvalues on the different datasets)
toplot <- tab[,c("sentinel","cohort", "graph_type", 
             "geuvadis_gene_gene","geo_gene_gene","cohort_gene_gene")]
colnames(toplot) <- c("sentinel", "cohort", "graph_type", 
                  "Geuvadis", "GEO", "LOLIPOP_KORA")

# plot the data
df <- melt(toplot,measure.vars=c(4,5,6), id.vars=c(1,2,3))
df$value <- -log10(df$value)
gp <- ggplot(df, aes(x=reorder(sentinel,value), y=value, fill=variable)) + 
  facet_grid(cohort+graph_type~variable) + sfm +
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x = element_text(vjust=1, angle = 90)) + 
  xlab("sentinel") +
  geom_hline(yintercept = -log10(0.05))

ggsave(plot=gp,
       file=fexpr,
       width=12, height=8)
