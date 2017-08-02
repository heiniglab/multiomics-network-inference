#'
#' The scripts loads a dataset and applies the bdgraph algorithm
#' to the data matrix, making use of predefined priors
#' 
#' @date 20170718
#' 
#' @author Johann Hawe
#'
library(BDgraph, lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4")

source("R/lib.R")
source("R/priors.R")

env <- new.env()
load("data/lolipop/rs9859077.data.RData", envir=env)
ggm.data <- with(env, sentinel)
id <- with(env, id)
data.matrix <- ggm.data$data

ranges <- ggm.data
ranges$data <- NULL

nodes <- colnames(data.matrix)
genes <- nodes[!grepl("^cg|^rs", nodes)]
# create edge matrix for all cpg.gene-cpg edges to be added
# to g.start
em <- matrix(ncol=2,nrow=0)
cpg.genes <- ranges$cpg.genes
cpg.genes <- cpg.genes[cpg.genes$SYMBOL %in% nodes]
cpgs <- (ranges$cpgs)
cpgs <- cpgs[names(cpgs) %in% nodes]
for(i in 1:length(cpgs)){
  cpg <- cpgs[i]
  g <- get.nearby.ranges(cpgs[i],cpg.genes)[[1]]
  
  for(row in 1:length(g)){
    em <- rbind(em, 
                c(names(cpg), 
                  g[row]$SYMBOL))  
  }
  
}

# create start graph
g.start <- create.g.start(nodes,
                          genes,
                          em)

# currently not implemented...
# priors <- get.link.priors(ranges)

ggm.fit <- bdgraph(data.matrix, 
        method="gcgm", 
        iter=100000, 
        burnin=20000,
        save.all=T, g.start = g.start)
save(file="results/ggm.fit.100k.full.RData", ggm.fit, ranges, data.matrix)
#traceplot(ggm.fit)
#plotcoda(ggm.fit)

g <- graphNEL.from.result(ggm.fit, 0.8, ranges)
plot.data <- plot.ggm(g=g, id=id)

# test a subset of the data matrix
samp <- sample(colnames(data.matrix), size = 20, replace=F)
if(!id %in% samp) samp <- c(samp,id)

data.test <- data.matrix[,samp]
priors.test <- priors[samp,samp]

ggm.test.fit <- bdgraph(data.test, 
        method="gcgm", 
        iter=3000, 
        save.all=T, 
        link.priors=priors.test)
graph <- graphNEL.from.result(ggm.test.fit, 0.9, ranges)
plot.data <- plot.ggm(graph,id=id, ranges)
save(file="results/rs9859077.test.fit.RData", 
     plot.data, 
     graph, 
     ranges,
     data.test,
     priors.test,
     ggm.test.fit)

traceplot(ggm.fit)
plotcoda(ggm.fit)
