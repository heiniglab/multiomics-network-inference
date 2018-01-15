# this script calls for each of the available sentinels (parsed from the data)
# the ggm modelling scripts.

# first load the custom BDgraph library plus others
library(BDgraph, lib.loc = "/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
library(pheatmap)
library(GenomicRanges)
library(igraph)
library(graph)

# get parameters
args <- commandArgs(trailingOnly = T)
sentinel <- args[1]
cores <- args[2]

plotdir <- "results/current/plots/"
outdir <- "results/current/fits/"

cat("Using sentinel", sentinel, "\n")
cat("Number of cores used:", cores, "\n")

#In this script we load KORA and LOLIPOP data for the `r sentinel` sentinel and 
#apply the bayesian GGM algorithm implemented in the BDgraph package to them, 
#using priors for the individual links possible in the underlying graph structure. 
#Priors have to be created beforehand and provided as RData files on the results folder.

# source the library scripts
source("R/lib.R")
source("R/priors.R")

# define the available cohorts
cohorts <- c("lolipop", "kora")

# load the sentinel data
load("results/current/data.processed.RData")
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

# set ggm parameters
iter=50000
burnin=40000

#We use `r iter` iterations with a burnin of `r burnin` as well as a total of
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
    
    # the file in which the model fit gets saved
    fit.out <- paste0(outdir, sentinel, ".", c, ".RData")
    if(!file.exists(fit.out)){
      ggm.fit <- bdgraph(gdata, 
                         method="gcgm", 
                         iter=iter, 
                         burnin=burnin,
                         save.all=T, g.start = gstart,
                         g.prior = gpriors, cores=cores)
      
      graph <- graph.from.fit(ggm.fit, ranges)
      if(length(nodes(graph))<1){
        warning("Resulting graph has no nodes.")
      } else {
        
        plot.data <- plot.ggm(g=graph, id=id, dot.out = paste0(outdir, sentinel, ".", c, ".dot"))
        save(file=fit.out, 
             ggm.fit, ranges, gdata, gstart, gpriors, graph, plot.data)
      }
    } else {
      print("Loading already available fit.")
      load(fit.out)
    }
    # check some plots
    pdf(paste0(plotdir, sentinel, ".", c, ".ggm.info.pdf"))  
    ggm.summary <- summary(ggm.fit)
    traceplot(ggm.fit)
    plotcoda(ggm.fit)
    dev.off()
    return(list(ggm.fit = ggm.fit, graph=graph, ranges=ranges, summary = ggm.summary))
  },
  #try-catch error
  error=function(m){
    f <- paste0(outdir, sentinel, ".out")
    cat("ERROR for sentinel", sentinel, "in cohort", c, "\n", 
        file=f)
    cat(paste0(m,collapse="\n"),file=f, append=T)
    cat(file="results/current/apply.ggm.failed.txt", append=T, sentinel, "\t", c)
    return(NULL)
  })
  fit
})

# check for error
if(!any(is.null(unlist(fits)))) {
  names(fits) <- cohorts
  # save fits extra, although this is a bit redundant
  save(file=paste0(outdir, sentinel, ".RData"), fits)
}
