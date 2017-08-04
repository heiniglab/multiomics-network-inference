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

# get the prior definitions
priors <- get.link.priors(ranges, nodes)

# create start graph
g.start <- get.g.start.from.priors(priors)

ggm.fit <- bdgraph(data.matrix, 
        method="gcgm", 
        iter=50000, 
        burnin=20000,
        multi.update = 1,
        save.all=T, g.start = g.start, g.prior = priors)

save(file="results/ggm.fit.100k.full.RData", 
     ggm.fit, ranges, priors, g.start, data.matrix)

g <- graphNEL.from.result(ggm.fit, 0.9, ranges)
plot.data <- plot.ggm(g=g, id=id, dot.out = "results/ggm.fit.100k.full.p09.dot")

q()

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
