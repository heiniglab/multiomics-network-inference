#' Script based on snakemake input/output vars to simulate
#' data for individual sentinels
#' 
#' @author Johann Hawe
#' 

# load needed packages
library(BDgraph, lib="~/epigenereg/packages/2017/R/3.4/")

# source needed scripts
source("R/priors.R")
source("R/lib.R")

# input file containing ranges/data of sentinel
ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]

# load data and utilize lolipop cohort
load(ifile)

for(s in "rs9859077") { #names(data)){
  ranges <- data[[s]]$lolipop$ranges
  nodes <- data[[s]]$lolipop$nodes
  ggm.data <- data[[s]]$lolipop$data
  
  cpgs <- with(ranges, names(cpgs))
  genes <- with(ranges, c(snp.genes$SYMBOL, cpg.genes$SYMBOL))
  if(!is.null(ranges$tfs)){
    genes <- c(genes, ranges$tfs$SYMBOL)
  }
  if(!is.null(ranges$spath)){
    genes <- c(genes, ranges$spath$SYMBOL)
  }
  genes <- unique(genes)

  # show some info
  print(paste0("Sentinel ID: ", s))
  print(paste0("Number of CpGs: ", length(cpgs)))
  print(paste0("Number of Genes: ", length(genes)))

  # get the priors (needed to create covariance matrix)
  priors <- get.link.priors(ranges, colnames(processed))
  
  # create a 'well-informed' covariance and precision matrix
  sigma <- priors
  # here we add some noise on the priors
  len <- dim(sigma)[1] * dim(sigma)[2]
  noise <- matrix(rnorm(len, 0, 1), dim(sigma)[1])
  sigma <- sigma + noise
  # diag should be all 1s, since this is cov(x,x)
  diag(sigma) <- 1
  K <- solve(sigma)
  # now we can simulate the data for the entities
  N <- nrow(ggm.data)
  p <- length(nodes)

  # define true graph
  # for now just use the start graph (i.e. where there was any
  # prior information available)
  G <- get.g.start.from.priors(priors)
  data.sim <- bdgraph.sim(p, graph=G, N, mean = 0)
  ggm.fit <- bdgraph(data.sim$data, g.start = G, method = "gcgm", 
                     iter = 5000, g.prior = priors, cores = 10)
}



# TODO create 'true' graph
g <- matrix()

