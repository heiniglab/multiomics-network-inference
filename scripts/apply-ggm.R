#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Prepare libraries and source scripts
# ------------------------------------------------------------------------------
# first load the custom BDgraph library plus others
library(BDgraph, lib.loc = "/home/icb/johann.hawe/R/x86_64-redhat-linux-gnu-library/3.4")
library(pheatmap)
library(GenomicRanges)
library(igraph)
library(graph)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# get snakemake parameters
# ------------------------------------------------------------------------------
franges <- snakemake@input[["ranges"]]
fdata <- snakemake@input[["data"]]
fstring <- snakemake@input[["string"]]
fpriors <- snakemake@input[["priors"]]
fcontext <- snakemake@input[["cpg_context"]]
fout <- snakemake@output$fit
fsummary_plot <- snakemake@output$summary_file
fgstart_plot <- snakemake@output$gstart_file

cores <- snakemake@threads
nriter <- snakemake@params$nriter
burnin <- snakemake@params$burnin

# ------------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------------
data <- readRDS(fdata)
priors <- readRDS(fpriors)
# filter for available data in priros
priors <- priors[colnames(data), colnames(data)]
ranges <- readRDS(franges)

# load string DB
STRING_DB <- readRDS(fstring)
  
# create the start graph for the GGM algorithm
gstart <- get_gstart_from_priors(priors)

# ------------------------------------------------------------------------------
# create the fits, learn our models
# ------------------------------------------------------------------------------
ggm_fit <- bdgraph(data, 
                   method="gcgm", 
                   iter=nriter, 
                   burnin=burnin,
                   save.all=T, 
                   g.start = gstart,
                   g.prior = priors, 
                   cores=cores)
ggm_fit_no_priors <- bdgraph(data, 
                             method="gcgm", 
                             iter=nriter, 
                             burnin=burnin,
                             save.all=T, 
                             g.start = gstart,
                             cores=cores)
      
# get the graphs from the ggm fits
graph <- graph_from_fit(ggm_fit, ranges, fcontext=fcontext)
graph_no_priors <- graph_from_fit(ggm_fit_no_priors, ranges, fcontext=fcontext)

# ------------------------------------------------------------------------------
# All done, finalize everything
# ------------------------------------------------------------------------------

# report some stats
print("Prior-graph:")
print(graph)
print("No-Prior-graph:")
print(graph_no_priors)

# create summary plots
pdf(fsummary_plot)  

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
# Save results
# ------------------------------------------------------------------------------
result <- listN(ggm_fit, ggm_fit_no_priors, graph, graph_no_priors, 
                gstart, ggm_summary, ggm_summary_no_priors)

saveRDS(file=fout, result)
