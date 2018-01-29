
sink(file=snakemake@log[[1]])

# load needed libraries
library(BDgraph, lib="~/epigenereg/packages/2017/R/3.4/")
library(GenomicRanges)

source("R/lib.R")
source("R/priors.R")

ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]
cores <- snakemake@threads
plotdir <- snakemake@params$plotdir

print(paste0("Processing file: ", ifile, 
             " using ", cores, " cores."))

# load data
load(ifile)

# set some parameters
iter <- 50000
burnin <- 10000

d <- data.sim$data
# sentinel name to generate plot filename
s <- colnames(d)[grepl("^rs", colnames(d))]

# get the priors
priors <- get.link.priors(ranges, colnames(d))
gstart <- get.g.start.from.priors(priors)

ggm_fit <- bdgraph(d, method = "gcgm",
                   iter = iter, burnin = burnin, 
                   g.prior = priors, g.start = gstart, 
                   save.all=T, cores=cores)

ggm_fit_no_priors <- bdgraph(d, method = "gcgm",
                             iter = iter, burnin = burnin, 
                             g.prior = priors, g.start = gstart, 
                             save.all=T, cores=cores)

# plot the bdgraph summary files
pdf(paste0(plotdir, s, "-simulation-ggmsummary.pdf"))
traceplot(ggm_fit)
plotcoda(ggm_fit)
traceplot(ggm_fit_no_priors)
plotcoda(ggm_fit_no_priors)
dev.off()

# get the result graph
ggm_graph <- graph.from.fit(ggm_fit, ranges)
ggm_graph_no_priors <- graph.from.fit(ggm_fit_no_priors, ranges)
print(paste0("Saving results to ", ofile))
save(file=ofile, ggm_fit, ranges, data.sim, ggm_graph, ggm_graph_no_priors)

sink()