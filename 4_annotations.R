## needed from previous script 2_edgeR_diff.R:
## a TopTags (edgeR) object listing genes tested for differential
## expression

## The Bioconductor mouse annotation libary
library(Mus.musculus)

## merge all this information in a data frame to combine mouse
## annotations with gene_id from Cufflinks get the information from
## this annotation dbi into data frames
.get.annot.frame <- function(){
    ens2geneID <- as.data.frame(org.Mm.egENSEMBL)
    name2SY <- as.data.frame(org.Mm.egSYMBOL2EG) # get 5 digit "symbol" to "short gene name" mapping
    SY2Lname<- as.data.frame(org.Mm.egGENENAME) # get 5 digit "symbol" to "long name name" mapping
    SY2GO <- as.data.frame(org.Mm.egGO) # get 5 digit "symbol" to gene ontology mapping
    annot.frame <- merge(name2SY, ens2geneID, by = "gene_id")
    annot.frame <- merge(annot.frame, SY2GO, by = "gene_id")
    annot.frame <- merge(annot.frame, SY2Lname, by = "gene_id")
    return(annot.frame)
}

annot.frame <- .get.annot.frame()

#####  objects created here for further use:
write.table(annot.frame, "output_data/annotation_data.csv", sep=",")



