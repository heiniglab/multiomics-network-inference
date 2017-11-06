get.gsc <- function(db, idtype) {
  ## idtype can by ENSEMBL, SYMBOL, ENTREZID, etc. see keytypes(db) for id types
  require(AnnotationDbi)
  require(GSEABase)
  
  frameData = select(db, keys=keys(db, idtype), keytype=idtype, columns=c("GO", "EVIDENCE", idtype))
  frameData = frameData[!is.na(frameData$EVIDENCE),]
  
  
  frame = GOFrame(frameData[,c("GO", "EVIDENCE", idtype)])
  allFrame = GOAllFrame(frame)
  
  gsc <- GeneSetCollection(allFrame, setType=GOCollection())
  return (gsc)
}

gsc.file = "data/current/geneset-collection.RData"
if (!file.exists(gsc.file)) {
  library(Homo.sapiens)
  gsc = get.gsc(Homo.sapiens, "SYMBOL")
  dir.create(dirname(gsc.file), recursive=T)
  save(gsc, file=gsc.file)
} else {
  load(gsc.file)
}


go.enrichment <- function(genes, universe, gsc, ontologies=c("MF", "BP", "CC")) {
  require(GSEABase)
  require(GOstats)

  go.tab = NULL
  for (ontology in ontologies) {
    params = GSEAGOHyperGParams(name="GO", ontology=ontology, geneIds=genes, universeGeneIds=universe, pvalueCutoff=1, testDirection="over", geneSetCollection=gsc, conditional=FALSE)
    hgt = hyperGTest(params)
    res = data.frame(ontology, summary(hgt, pvalue=1))
    res = data.frame(res, q=p.adjust(res[,"Pvalue"]))
    colnames(res) = gsub("^GO.*ID$", "GOID", colnames(res))
    go.tab = rbind(go.tab, res)
  }
  return(go.tab)
}
