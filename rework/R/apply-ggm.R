#!/usr/bin/env Rscript

# this script calls for each of the available sentinels (parsed from the data)
# the ggm modelling scripts.

log <- snakemake@log[[1]]
sink(log)

# first load the custom BDgraph library plus others
library(BDgraph, lib.loc = "/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
library(pheatmap)
library(GenomicRanges)
library(igraph)
library(graph)

# get parameters
sentinel <- snakemake@wildcards$sentinel
cores <- snakemake@threads

# whether no priors should used in GGM
# we always use priors for gstart generation
#nopriors <- snakemake@params

if(is.null(sentinel)) {
  stop("No sentinel id provided.")
}

outfile <- snakemake@output[[1]]
plotdir <- snakemake@params$plotdir
nriter <- snakemake@params$nriter
burnin <- snakemake@params$burnin

cat("Using sentinel", sentinel, "\n")
cat("Number of cores used:", cores, "\n")

#In this script we load KORA and LOLIPOP data for the `r sentinel` sentinel and 
#apply the bayesian GGM algorithm implemented in the BDgraph package to them, 
#using priors for the individual links possible in the underlying graph structure. 
#Priors have to be created beforehand and provided as RData files on the results folder.

# source the library scripts
source("R/lib.R")
source("R/priors.R")
source("R/bdgraph-supplement.R")

# define the available cohorts
cohorts <- c("lolipop", "kora")

# load the sentinel data
load(snakemake@input[[1]])
data <- data[[sentinel]]

## Prior calculation
#Now we get the priors for all possible links within the data matrix. The priors
#need to have been created beforehand and are simply loaded in the called function.

priors <- lapply(cohorts, function(c){
  cohort <- data[[c]]

  # get the prior definitions and plot the prior heatmap
  priors <- get.link.priors(cohort$ranges, cohort$nodes)
  cat("Prior min-value: ", min(priors), "\n")
  print(dim(priors))
  pheatmap(priors, cex=0.7, main=paste0(c, " priors"),
           filename=paste0(plotdir, sentinel, ".", c, ".priors.pdf"),
           cex=0.7)
   
  return(priors)
})
names(priors) <- cohorts

## GGM fit
#In the next step we finally fit the GGM to our data, using the priors during 
#the process. First a 'start-graph' is created, which serves as the entry point
#in the calculations for the algorithms. Afterwards we can start running the 
#bdgraph algorithm in order to get our partial correlations based on data and prior knowledge. 

gstarts <- lapply(cohorts, function(c){
  # create start graph
  g.start <- get.g.start.from.priors(priors[[c]])
  pheatmap(g.start, cex=0.7, main=paste0(c, " start graph"),
           filename=paste0(plotdir, sentinel, ".", c, ".gstart.pdf"),
           cex=0.7)
  
  return(g.start)
})
names(gstarts) <- cohorts

#We use `r nriter` iterations with a burnin of `r burnin` as well as a total of
#`r cores` cores. We then extract the graph from the ggm fit and create a dot file
#representing the graph. If the graph is less then 500 edges it will be plotted graphically
#as well.
#In addition, we create some diagnostic plots to check the convergence of the 
#algorithm.


fits <- lapply(cohorts, function(c) {
  fit <- tryCatch({
    gdata <- data[[c]]$data
    id <- data[[c]]$id
    ranges <- data[[c]]$ranges
    gstart <- gstarts[[c]]
    gpriors <- priors[[c]]

    # create the fits
    ggm_fit <- bdgraph(gdata, 
                         method="gcgm", 
                         iter=nriter, 
                         burnin=burnin,
                         save.all=T, g.start = gstart,
                         g.prior = gpriors, cores=cores)
    ggm_fit_no_priors <- bdgraph(gdata, 
                                 method="gcgm", 
                                 iter=nriter, 
                                 burnin=burnin,
                                 save.all=T, g.start = gstart,
                                 cores=cores)
      
    graph <- graph.from.fit(ggm_fit, ranges)
    graph_no_priors <- graph.from.fit(ggm_fit_no_priors, ranges)
    if(length(nodes(graph))<1){
      warning("Resulting graph has no nodes.")
    } else {        
      plot.data <- plot.ggm(g=graph, id=id, dot.out = gsub("\\..+", ".dot", outfile))
    }

    # check some plots
    pdf(paste0(plotdir, sentinel, ".", c, ".ggm.info.pdf"))  
    ggm_summary <- summary(ggm_fit)
    traceplot(ggm_fit)
    plotcoda(ggm_fit)
    ggm_summary_no_priors <- summary(ggm_fit_no_priors)
    traceplot(ggm_fit_no_priors)
    plotcoda(ggm_fit_no_priors)
    dev.off()
    return(list(ggm_fit, ggm_fit_no_priors, graph, graph_no_priors, ranges, 
                gstart, gpriors, gdata, ggm_summary, ggm_summar_no_priors))
  },
  #try-catch error
  error=function(m){
    message(paste0("ERROR for sentinel ", sentinel, " in cohort ", c, "\n"))
    message(paste0(m,collapse="\n"))
    return(NULL)
  })
  fit
})

# check for error
if(!any(is.null(unlist(fits)))) {
  names(fits) <- cohorts
  save(file=outfile, fits)
}
sink()
