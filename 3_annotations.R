## needed from previous script 2_edgeR_diff.R:
## a TopTags (edgeR) object listing genes tested for differential
## expression

if(!exists("ALL.top")){
    source("2_edgeR_diff.R")
}

## The Bioconductor mouse annotation libary
library(Mus.musculus)

get.cuff2name <- function(path="data/merged_mm9_Ef.gtf"){
    gene.gtf.info <- gffRead(path) 
    gff.attr.list <- strsplit(gene.gtf.info$attributes, " ")
    cuff_ids <- unlist(lapply(gff.attr.list, function (x) x[[4]]))
    gene_names <- unlist(lapply(gff.attr.list, function (x) x[[6]]))
    cuff_ids <- gsub("\"(.*)\";$", "\\1", cuff_ids)
    gene_names <- gsub("\"(.*)\";$", "\\1", gene_names)
    ## discard the cufflinks re-annotation information
    gene_names <- gsub("_\\d$", "", gene_names) 
    cuff2name <-as.data.frame(cbind(cuff_ids, gene_names))
    return(unique(cuff2name))
}

cuff2name <- get.cuff2name()

## merge all this information in a data frame to combine mouse
## annotations with gene_id from Cufflinks get the information from
## this annotation dbi into data frames
.get.annot.frame <- function(){
    name2SY <- as.data.frame(org.Mm.egSYMBOL2EG) # get 5 digit "symbol" to "short gene name" mapping
    SY2Lname<- as.data.frame(org.Mm.egGENENAME) # get 5 digit "symbol" to "long name name" mapping
    SY2GO <- as.data.frame(org.Mm.egGO) # get 5 digit "symbol" to gene ontology mapping
    annot.frame <- merge(name2SY, cuff2name, by.x = "symbol", by.y = "gene_names")
    annot.frame <- merge(annot.frame, SY2GO, by = "gene_id")
    annot.frame <- merge(annot.frame, SY2Lname, by = "gene_id")
    return(annot.frame)
}

annot.frame <- .get.annot.frame()

## cuff2GO contains the mapping of cufflink ids to GO terms
cuff2GO.Mm <- by(annot.frame, annot.frame$cuff_ids,
                 function(x) as.character(x$go_id))
## define the universe of tested genes
cuff2GO.Mm <- cuff2GO.Mm[names(cuff2GO.Mm)%in%rownames(ALL.top.Mm)]

## overwrite the cuff2name creator function to turn it into a subsetter function
get.cuff2name <- function(cuff_ids){
    cuffsubset <- cuff2name[cuff2name$cuff_ids%in%cuff_ids , "gene_names"]
    return(as.character(cuffsubset))
}

go_Ef <- read.csv("data/12864_2014_6777_GO-annot.csv")

gene2GO.Ef <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))

## subset for only those genes known that have been tested to create
## the suitable universe to search enrichment against


## objects created here for further use:

## 1. cuff2name - a data frame linking gene ids
## 2. annot.frame - a data frame containing coprehensive annotation information
## 3. cuff2GO.Mm - a list linking cufflinks ids to GO terms
## 4. get.cuff2name - a function for extracting proper (mouse or Efalciformis) gene ids for a cufflink id
