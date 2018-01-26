
sink(file=snakemake@log[[1]])

# load needed libraries
library(BDgraph, lib="~/epigenereg/packages/2017/R/3.4/")
library(GenomicRanges)

source("R/lib.R")
source("R/priors.R")

ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]
cores <- snakemake@threads

print(paste0("Processing file: ", ifile, 
             " using ", cores, " cores."))

# load data
load(ifile)

# set some parameters
iter <- 10000
burnin <- 5000

d <- data.sim$data

# get the priors
priors <- get.link.priors(ranges, colnames(d))
gstart <- get.g.start.from.priors(priors)

ggm_fit <- bdgraph(d, method = "gcgm",
                   iter = iter, burnin = burnin, 
                   g.prior = priors, g.start = gstart, 
                   save.all=T, cores=)

# get the result graph
ggm_graph <- graph.from.fit(ggm_fit, ranges)

print(paste0("Saving results to ", ofile))
save(file=ofile, ggm_fit, ranges, data.sim)

sink()