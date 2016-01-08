# Annotations Eimeria falciformis

#name2SY <- as.list(org.Mm.egSYMBOL2EG) # get 5 digit symbol from XLOC
#SY2Lname<- as.list(org.Mm.egGENENAME) # use 5 digit symbol to get long name for gene
#SY2GO <- as.list(org.Mm.egGO) # get gene ontology annotations
#SY2GO <- SY2GO[!is.na(SY2GO)]

# csv.files taken from Heitlinger et al. on Ef genome
#gene.gtf.info <- gffRead("/data/Eimeria_Totta/reference_genomes/merged_mm9_Ef.gtf") 
# gffRead() not in available libraries - could not find it - hand-built?
#gff.attr.list <- strsplit(gene.gtf.info$attributes, " ")

#gene_ids <- unlist(lapply(gff.attr.list, function (x) x[[4]]))
#gene_names <- unlist(lapply(gff.attr.list, function (x) x[[6]]))

#gene_ids <- gsub("\"(.*)\";$", "\\1", gene_ids)
#gene_names <- gsub("\"(.*)\";$", "\\1", gene_names)

#id2name <-as.data.frame(cbind(gene_ids, gene_names))

#id2name <- unique(id2name)


############ below line does not worki#################
#flat.GO <- lapply(SY2GO, function (x) do.call("rbind", x))

#GO.df <- as.data.frame(do.call("rbind", flat.GO))
#GO.df$gene_id <- rep(names(flat.GO), times = unlist(lapply(flat.GO, nrow)))

gene2GO <- by(GO.df, GO.df$gene_id, function(x) as.character(x$GOID))

## subset for only those genes known that have been tested to create
## the suitable universe to search enrichment against

all.genes <- unique(select.from.stats.results(stat_results_strain, mouse.bg,
                                              qval = 2))

all.genes <- get.annotation.for.xloc(all.genes)[[2]]
gene2GO <- gene2GO[names(gene2GO)%in%all.genes]

