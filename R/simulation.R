library(graph)
library(BDgraph, lib="~/R/x86_64-redhat-linux-gnu-library/3.4/")
library(mvtnorm)
library(pheatmap)
library(RColorBrewer)

# number of iterations for bdgraph
bd.iter <- 10000;
# number of simulations
sim <- 1000
# output dir for plot
ODIR <- "results/simulation"
dir.create(ODIR)

cat("Started processing at ", date(), "\n")

#' Update counts
#'
#' Specifically updates the counts for an edgelist yielded during the testing of
#' the bdgraph implementation.
#'
#' @param n The name in the 'counts' variable for which to update its counts
#' @param elist The list of edges obtained from the respective graphNEL object
#'
update.counts <- function(n, elist) {
  temp <- counts[[n]]
  # check which edges are in the graph and increase our overall counts
  if("x1" %in% names(elist)) {
    if("x2" %in% elist[["x1"]]) {
      temp["x1x2"] <- temp["x1x2"] + 1;
    }
    if("x3" %in% elist[["x1"]]) {
      temp["x1x3"] <- temp["x1x3"] + 1;
    }
    if("x4" %in% elist[["x1"]]) {
      temp["x1x4"] <- temp["x1x4"] + 1;
    }
  }
  if("x2" %in% names(elist)){
    if("x3" %in% elist[["x2"]]) {
      temp["x2x3"] <- temp["x2x3"] + 1;
    }
    if("x4" %in% elist[["x2"]]) {
      temp["x2x4"] <- temp["x2x4"] + 1;
    }
  }
  if("x3" %in% names(elist)){
    if("x4" %in% elist[["x3"]]) {
      temp["x3x4"] <- temp["x3x4"] + 1;
    }
  }
  counts[[n]] <<- temp;

}

# counts how often a certain edge was observed
x1x2 <- 0
x1x3 <- 0
x1x4 <- 0
x2x3 <- 0
x2x4 <- 0
x3x4 <- 0

counts <- list(e1=c(x1x2=x1x2,x1x3=x1x3,x1x4=x1x4,x2x3=x2x3,x2x4=x2x4,x3x4=x3x4),
               e2=c(x1x2=x1x2,x1x3=x1x3,x1x4=x1x4,x2x3=x2x3,x2x4=x2x4,x3x4=x3x4));

# number of nodes
p <- 4
# number of samples
n <- 400
# covariance between samples to be used
covs <- c(0.1,0.3,0.6)

for(co in covs){
  # reset counts between iterations
  counts <- list(e1=c(x1x2=x1x2,x1x3=x1x3,x1x4=x1x4,x2x3=x2x3,x2x4=x2x4,x3x4=x3x4),
               e2=c(x1x2=x1x2,x1x3=x1x3,x1x4=x1x4,x2x3=x2x3,x2x4=x2x4,x3x4=x3x4))
  cat("Cov:", co, "\n")
  # main loop for simulation
  for(i in 1:sim) {
    cat("sim",i,"\n")
    dat <- rmvnorm(n, sigma = matrix(c(1, rep(co,3),
                                       co,1, rep(co,2),
                                       rep(co,2),1,co,
                                       rep(co,3),1), ncol=p));
  
    colnames(dat) <- paste0("x", 1:p)
    # create priors
    priors <- matrix(0.0000001,ncol=p,nrow=p)
    colnames(priors) <- rownames(priors) <- paste0("x", 1:p)
  
    priors[2,3] <- priors[3,2] <- 0.9
    priors[3,4] <- priors[4,3] <- 0.5
  
    bdtest <- bdgraph(dat, iter=bd.iter, method="ggm", print=floor(bd.iter/2),
                      save.all=T, g.prior=priors, cores=20)
  
    g.adj1 <- summary(bdtest, vis=F)$selected_g
    g.adj2 <- BDgraph::select(bdtest, cut = 0.9)
  
    graph1 <- as_graphnel(graph.adjacency(g.adj1, mode="undirected", diag=F));
    graph2 <- as_graphnel(graph.adjacency(g.adj2, mode="undirected", diag=F));
  
    e1 <- graph::edges(graph1)
    e2 <- graph::edges(graph2)
    update.counts("e1",e1)
    update.counts("e2",e2)
  
  }
  
  
  # make simple barplots displaying the results
  pdf(paste0(ODIR, "/", sim,"simulations_", p, "p_", co, "cov.pdf"))
  m <- paste0("Number of edges between nodes aggregated over ", sim," simulations.")
  cols <- brewer.pal(p*(p-1)/2, "Dark2")
  barplot(counts[["e1"]], main=m, col=cols)
  legend("topleft", legend=names(counts[["e1"]]), fill=cols)
  barplot(counts[["e2"]], main=m, col=cols)
  legend("topleft", legend=names(counts[["e2"]]), fill=cols)
  dev.off()

}
cat("Finished processing at ", date(), "\n")

sessionInfo()

