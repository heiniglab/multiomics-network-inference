#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
#' Apply graphical model inference on data collected for a specific sentinel.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Prepare libraries and source scripts.")
# ------------------------------------------------------------------------------
# first load the custom BDgraph library plus others
library(BDgraph,
        lib.loc = "/home/icb/johann.hawe/R/x86_64-redhat-linux-gnu-library/3.4")
library(pheatmap)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(graph))
library(GeneNet)
source("scripts/lib.R")
source("scripts/reg_net.R")

# ------------------------------------------------------------------------------
print("Get snakemake parameters.")
# ------------------------------------------------------------------------------

#input
franges <- snakemake@input[["ranges"]]
fdata <- snakemake@input[["data"]]
fppi_db <- snakemake@input[["ppi_db"]]
fpriors <- snakemake@input[["priors"]]
fcontext <- snakemake@input[["cpg_context"]]

# output
fout <- snakemake@output$fit
fsummary_plot <- snakemake@output$summary_file
fgstart_plot <- snakemake@output$gstart_file

# params
cores <- snakemake@threads
nriter <- snakemake@params$nriter
burnin <- snakemake@params$burnin

# ------------------------------------------------------------------------------
print("Load and prepare data.")
# ------------------------------------------------------------------------------
data <- readRDS(fdata)
priors <- readRDS(fpriors)

# filter for available data in priros
priors <- priors[colnames(data), colnames(data)]
ranges <- readRDS(franges)

# load PPI DB
ppi_db <- readRDS(fppi_db)

# create the start graph for the GGM algorithm
gstart <- get_gstart_from_priors(priors)

# ------------------------------------------------------------------------------
# for catching the genenet summary plots, we open the pdf connection here
# ------------------------------------------------------------------------------
pdf(fsummary_plot)

# ------------------------------------------------------------------------------
print("Infer regulatory networks.")
# ------------------------------------------------------------------------------
print("Fitting model using priors.")
bdgraph <- reg_net(d, priors, "bdgraph", threads=threads)

#print("Fitting model with priors without start graph")
#bdgraph_empty <- reg_net(d, priors, "bdgraph",
#                         use_gstart = F, threads=threads)

#print("Fitting model without priors with start graph")
#bdgraph_no_priors <- reg_net(d, NULL, "bdgraph",
#                             gstart = gstart, threads=threads)

print("Fitting model without priors using empty start graph.")
bdgraph_no_priors_empty <- reg_net(d, NULL, "bdgraph",
                                   use_gstart = F, threads=threads)

print("Fitting model using iRafNet.")
irafnet <- reg_net(d, priors, "irafnet", threads=threads)

print("Fitting model using GeneNet.")
genenet <- reg_net(d, priors, "genenet", threads=threads)

# ------------------------------------------------------------------------------
print("Add custom annotations for the graphs.")
# ------------------------------------------------------------------------------
bdgraph$graph <- annotate.graph(bdgraph$graph, ranges, ppi_db, fcontext)
bdgraph_no_priors_empty$graph <- annotate.graph(bdgraph_no_priors_empty$graph,
                                                ranges, ppi_db, fcontext)
irafnet$graph <- annotate.graph(irafnet$graph, ranges, ppi_db, fcontext)
genenet$graph <- annotate.graph(genenet$graph, ranges, ppi_db, fcontext)

# ------------------------------------------------------------------------------
print("Create result list.")
# ------------------------------------------------------------------------------
result <- list(ggm_fit = bdgraph$fit,
               ggm_fit_no_priors_empty = bdgraph_no_priors_empty$fit,
               irn_fit = irafnet$fit,
               genenet_fit = genenet$fit,
               ggm_graph = bdgraph$graph,
               ggm_graph_no_priors_empty = bdgraph_no_priors_empty$graph,
               irn_graph = irafnet$graph,
               genenet_graph = genenet$graph)

# ------------------------------------------------------------------------------
print("Done with model fitting. Plotting some information.")
# ------------------------------------------------------------------------------

# plot convergence info
ggm_summary <- summary(ggm_fit)
traceplot(ggm_fit)
plotcoda(ggm_fit)
ggm_summary_no_priors <- summary(ggm_fit_no_priors)
traceplot(ggm_fit_no_priors)
plotcoda(ggm_fit_no_priors)

# plot the start graph as a matrix
pheatmap(gstart, cex=0.7, main="start graph",
         filename=fgstart_plot,
         cex=0.7)

dev.off()

# ------------------------------------------------------------------------------
print("Plotting done. Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, result)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
