#' Performs validation of already fitted GGMs
#' 
#' @author Johann Hawe
#' 
#' @version 20170720
#'

source("R/validation.R")

# prepare data
load("results/rs9859077.test.fit.RData")
data <- data.test
snp <- names(ranges$sentinel.range)
data[,snp] <- as.integer(as.character(data[,snp]))
cpgs <- names(ranges$cpgs)
genes <- ranges$snp.genes$SYMBOL
genes <- genes[genes %in% colnames(data)]

gn <- nodes(graph)
selected.genes <- gn[unlist(nodeData(graph,gn,"snp.gene"))]

par(mfrow=c(4,round(length(genes)/4+0.5)), mar=rep(4.5,4))
pdf("senp7.mediation.pdf")
m <- mediation(data, snp, genes, cpgs)
dev.off()

# snp~gene link validation
vsnp <- validate.snps(ggm.fit, ranges)

# cpg~gene link validation
vcpg <- validate.cpgs(ggm.fit, ranges)

# gene~gene link validation
vgene <- validate.genes(ggm.fit, ranges)


