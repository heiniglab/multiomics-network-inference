
# load needed libraries
library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
library(GenomicRanges)

source("scripts/lib.R")
source("scripts/bdgraph-supplement.R")

fdata <- snakemake@input[[1]]
fout <- snakemake@output[[1]]
fplot <- snakemake@output[[2]]
cores <- snakemake@threads
nriter <- as.numeric(snakemake@params$nriter)
burnin <- as.numeric(snakemake@params$burnin)

# load data
# contains: simulated_data, ranges, nodes, ggm.data
load(fdata)

# we generated several graphs, which we all want to process now
result <- lapply(names(simulated_data), function(n) {
  sim <- simulated_data[[n]]
  # sentinel name to generate plot filename
  s <- sim$snp
  d <- sim$data.sim
  
  # get the priors and create start graph
  priors <- sim$priors
  # Start from observed/hidden graph state?
  gstart <- get.g.start.from.priors(priors)
  
  ggm_fit <- bdgraph(d, method = "gcgm",
                     iter = nriter, burnin = burnin, 
                     g.prior = priors, g.start = gstart, 
                     save.all=T, cores=cores)
  
  ggm_fit_no_priors <- bdgraph(d, method = "gcgm",
                               iter = nriter, burnin = burnin, 
                               g.start = gstart, 
                               save.all=T, cores=cores)
  
  # plot the bdgraph summary files
  pdf(fplot)
  temp <- summary(ggm_fit)
  plotcoda(ggm_fit)
  temp <- summary(ggm_fit_no_priors)
  plotcoda(ggm_fit_no_priors)
  dev.off()
  
  # get the result graph
  ggm_graph <- graph.from.fit(ggm_fit, ranges, annotate=F)
  ggm_graph_no_priors <- graph.from.fit(ggm_fit_no_priors, ranges, annotate=F)
 
  # new entry in our data collection
  sim$fits <- listN(ggm_fit, ggm_fit_no_priors, 
                   ggm_graph, ggm_graph_no_priors)
  sim
})

names(result) <- names(simulated_data)

save(file=fout, result)
