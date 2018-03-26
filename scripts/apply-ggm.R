#!/usr/bin/env Rscript

# this script calls for each of the available sentinels 

# first load the custom BDgraph library plus others
library(BDgraph, lib.loc = "/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
library(pheatmap)
library(GenomicRanges)
library(igraph)
library(graph)

# get parameters
franges <- snakemake@input[["ranges"]]
fdata <- snakemake@input[["data"]]
fpriors <- snakemake@input[["priors"]]
cores <- as.numeric(snakemake@threads)
fout <- snakemake@output[[1]]
fplot <- snakemake@params$plot_file
nriter <- as.numeric(snakemake@params$nriter)
burnin <- as.numeric(snakemake@params$burnin)

#In this script we load KORA and LOLIPOP data for the `r sentinel` sentinel and 
#apply the bayesian GGM algorithm implemented in the BDgraph package to them, 
#using priors for the individual links possible in the underlying graph structure. 
#Priors have to be created beforehand and provided as RData files on the results folder.

# source the library scripts
source("scripts/lib.R")
source("scripts/priors.R")
source("scripts/bdgraph-supplement.R")

data <- readRDS(fdata)
priors <- readRDS(fpriors)

## GGM fit
#In the next step we finally fit the GGM to our data, using the priors during 
#the process. First a 'start-graph' is created, which serves as the entry point
#in the calculations for the algorithms. Afterwards we can start running the 
#bdgraph algorithm in order to get our partial correlations based on data and prior knowledge. 

gstart <- get.g.start.from.priors(priors)

#We use `r nriter` iterations with a burnin of `r burnin` as well as a total of
#`r cores` cores. We then extract the graph from the ggm fit and create a dot file
#representing the graph. If the graph is less then 500 edges it will be plotted graphically
#as well.
#In addition, we create some diagnostic plots to check the convergence of the 
#algorithm.

# create the fits
ggm_fit <- bdgraph(data, 
                   method="gcgm", 
                   iter=nriter, 
                   burnin=burnin,
                   save.all=T, g.start = gstart,
                   g.prior = priors, cores=cores)
ggm_fit_no_priors <- bdgraph(data, 
                             method="gcgm", 
                             iter=nriter, 
                             burnin=burnin,
                             save.all=T, g.start = gstart,
                             cores=cores)
      
# get the graphs from the ggm fits
graph <- graph.from.fit(ggm_fit, ranges)
graph_no_priors <- graph.from.fit(ggm_fit_no_priors, ranges)

# report some stats
print("Prior-graph:")
print(graph)
print("No-Prior-graph:")
print(graph_no_priors)

# create summary plots
pdf(fplot)  

# plot convergence info
ggm_summary <- summary(ggm_fit)
traceplot(ggm_fit)
plotcoda(ggm_fit)
ggm_summary_no_priors <- summary(ggm_fit_no_priors)
traceplot(ggm_fit_no_priors)
plotcoda(ggm_fit_no_priors)
# plot the start graph as a matrix
pheatmap(g_start, cex=0.7, main=paste0(c, " start graph"),
         filename=paste0(plotdir, sentinel, ".", c, ".gstart.pdf"),
         cex=0.7)

dev.off()

result <- listN(ggm_fit, ggm_fit_no_priors, graph, graph_no_priors, 
                gstart, ggm_summary, ggm_summary_no_priors))

saveRDS(file=fout, result)
