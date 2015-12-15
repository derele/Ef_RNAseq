# Bioconductor mouse annotation libaries mouse

library(Mus.musculus)

##
## get the stuff from annotation dbi
name2SY <- as.list(org.Mm.egSYMBOL2EG) # get 5 digit symbol from XLOC
SY2Lname<- as.list(org.Mm.egGENENAME) # use 5 digit symbol to get long name for gene
SY2GO <- as.list(org.Mm.egGO) # get gene ontology annotations


## Combine mouse annotations with gene_id from Cufflinks output format (XLOC_)
gene.gtf.info <- gffRead("/data/Eimeria_Totta/reference_genomes/merged_mm9_Ef.gtf") 
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





