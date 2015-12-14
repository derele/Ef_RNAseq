# Bioconductor mouse annotation libaries mouse

library(Mus.musculus)
library(topGO)

##
if(!exists("chal_genes")) {
  source("2b_ballgown_diff_EH.R")
} 

## get the stuff from annotation dbi
name2SY <- as.list(org.Mm.egSYMBOL2EG) # get 5 digit symbol from XLOC
SY2Lname<- as.list(org.Mm.egGENENAME) # use 5 digit symbol to get long name for gene
SY2GO <- as.list(org.Mm.egGO) # get gene ontology annotations


## Combine mouse annotations with gene_id from Cufflinks output format (XLOC_)
gene.gtf.info <- gffRead("../reference_genomes/merged_mm9_Ef.gtf") # path assuming location in 'Rfolder'
gff.attr.list <- strsplit(gene.gtf.info$attributes, " ")

gene_ids <- unlist(lapply(gff.attr.list, function (x) x[[4]]))
gene_names <- unlist(lapply(gff.attr.list, function (x) x[[6]]))

gene_ids <- gsub("\"(.*)\";$", "\\1", gene_ids)
gene_names <- gsub("\"(.*)\";$", "\\1", gene_names)

id2name <-as.data.frame(cbind(gene_ids, gene_names))

id2name <- unique(id2name)


## Use this function to get annotations for any XLOC mouse gene-id
get.annotation.for.xloc <- function(xloc){ 
  sig.names <- as.character(id2name[id2name$gene_ids%in%xloc, "gene_names"])
  sig.SY <- unlist(name2SY[sig.names])
  sig.Lnames  <- unlist(SY2Lname[sig.SY])
  return(list(sig.names, sig.SY, sig.Lnames))
}


SY2GO <- SY2GO[!is.na(SY2GO)]
flat.GO <- lapply(SY2GO, function (x) do.call("rbind", x))


GO.df <- as.data.frame(do.call("rbind", flat.GO))
GO.df$gene_id <- rep(names(flat.GO), times = unlist(lapply(flat.GO, nrow)))


gene2GO <- by(GO.df, GO.df$gene_id, function(x) as.character(x$GOID))

## subset for only those genes known that have been tested to create
## the suitable universe to search enrichment against

all.genes <- unique(select.from.stats.results(stat_results_strain, mouse.bg,
                                              qval = 2))

all.genes <- get.annotation.for.xloc(all.genes)[[2]]
gene2GO <- gene2GO[names(gene2GO)%in%all.genes]


TOGO.all.onto <- function (ontology, allgenes, gene.set, annot) {
  g <- factor(as.integer( allgenes%in%gene.set ))
  names(g) <- allgenes
  toGO <-  new("topGOdata", ontology = ontology, allGenes = g, annot = annFUN.gene2GO,
               gene2GO = annot)
  resultFis <- runTest(toGO, algorithm = "classic", statistic = "fisher")
  list(toGO, resultFis) ## returns a list first data then result
}


gene.table.topGO <- function(TOGO.list, pval=0.01){
    all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
    ##    all$fdr <- p.adjust(all$result1, method="BH")
    names(all)[names(all)%in%"result1"] <- "p.value"
    return(all[all$p.value<pval,])
}



chalGids <- get.annotation.for.xloc(chal_genes)[[2]]

MF.chal <- TOGO.all.onto("MF", names(gene2GO), 
                         chalGids, gene2GO)
gene.table.topGO(MF.chal)

BP.chal <- TOGO.all.onto("BP", names(gene2GO), 
                         chalGids , gene2GO)
head(gene.table.topGO(BP.chal), n=20)



chalGidsI <- get.annotation.for.xloc(chal_genesI)[[2]]

MF.chalI <- TOGO.all.onto("MF", names(gene2GO), 
                         chalGidsI, gene2GO)
gene.table.topGO(MF.chalI)

BP.chalI <- TOGO.all.onto("BP", names(gene2GO), 
                          chalGidsI , gene2GO)
head(gene.table.topGO(BP.chalI, 1), n=20)



strainGids5 <- get.annotation.for.xloc(strain_genes5)[[2]]

MF.strain5 <- TOGO.all.onto("MF", names(gene2GO), 
                            strainGids5, gene2GO)
gene.table.topGO(MF.strain5)

BP.strain5 <- TOGO.all.onto("BP", names(gene2GO), 
                            strainGids5, gene2GO)
head(gene.table.topGO(BP.strain5), n=20)


strainGidsI <- get.annotation.for.xloc(strain_genesI)[[2]]

MF.strainI <- TOGO.all.onto("MF", names(gene2GO), 
                            strainGidsI, gene2GO)
gene.table.topGO(MF.strainI)

BP.strainI <- TOGO.all.onto("BP", names(gene2GO), 
                            strainGidsI, gene2GO)
head(gene.table.topGO(BP.strainI), n=20)






