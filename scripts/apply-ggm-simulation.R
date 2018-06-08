# ------------------------------------------------------------------------------
# load needed libraries and source scripts
# ------------------------------------------------------------------------------
library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
library(GenomicRanges)

source("scripts/lib.R")
source("scripts/bdgraph-supplement.R")

# ------------------------------------------------------------------------------
# Get snakemake input and load data
# ------------------------------------------------------------------------------
fdata <- snakemake@input[[1]]
fout <- snakemake@output[[1]]
cores <- snakemake@threads
nriter <- as.numeric(snakemake@params$nriter)
burnin <- as.numeric(snakemake@params$burnin)
iter <- as.numeric(snakemake@params$iteration)

# contains: simulations, ranges, priors, nodes, data, run
load(fdata)

# ------------------------------------------------------------------------------
# we generated several graphs, which we all process now (i.e. calculate ggms)
# ------------------------------------------------------------------------------

# plot file for the bdgraph summary
#pdf(fplot)

# apply over the different runs/iterations
simulated_data <- simulations[[iter]]
  
# for this run, apply over all simulated graphs (different randomization
# degrees)
result <- lapply(names(simulated_data), function(n) {
  # --------------------------------------------------------------------------
  # Get data
  # --------------------------------------------------------------------------
  sim <- simulated_data[[n]]
  # sentinel name to generate plot filename
  s <- sim$snp
  d <- sim$data.sim
  
  # in case of rbinom simulation, we get adjusted priors
  if(grepl("_rbinom", n)) {
    priors <- sim$priors
  }
  
  # Start from observed/hidden graph state?
  gstart <- get_gstart_from_priors(priors)
  
  # --------------------------------------------------------------------------
  # Get model fits
  # --------------------------------------------------------------------------
  print("Fitting model using priors.")
  ggm_fit <- bdgraph(d, method = "gcgm",
                     iter = nriter, burnin = burnin, 
                     g.prior = priors, g.start = gstart, 
                     save.all=T, cores=cores)
  print("Fitting model without priors")
  ggm_fit_no_priors <- bdgraph(d, method = "gcgm",
                               iter = nriter, burnin = burnin, 
                               g.start = gstart, 
                               save.all=T, cores=cores)
    
  # --------------------------------------------------------------------------
  # Plot summaries
  # --------------------------------------------------------------------------
  #temp <- summary(ggm_fit)
  #plotcoda(ggm_fit, main=paste0("CODA:", n))
  #temp <- summary(ggm_fit_no_priors)
  #plotcoda(ggm_fit_no_priors, main=paste0("CODA:", n, " no priors"))
    
  # --------------------------------------------------------------------------
  # Extract graphs and create result lists
  # --------------------------------------------------------------------------
  # get the result graph
  ggm_graph <- graph_from_fit(ggm_fit, ranges, annotate=F)
  ggm_graph_no_priors <- graph_from_fit(ggm_fit_no_priors, ranges, annotate=F)
   
  # new entry in our data collection
  sim$fits <- listN(ggm_fit, ggm_fit_no_priors, 
                    ggm_graph, ggm_graph_no_priors)
  sim  
})
names(result) <- names(simulated_data)
#dev.off()

print("Saving results.")
save(file=fout, result, priors, runs)
