## needed from previous script 2_edgeR_diff.R:
## a TopTags (edgeR) object listing genes tested for differential
## expression

if(!exists("Ef.1st.pass.model")){ #Totta changes Mm to Ef here
    source("2_edgeR_diff.R")
}

# Check checkpoint package for package control.
if(!require(Mus.musculus)) biocLite("Mus.musculus") # Imports package if user does not have it

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

## ensembl2GO contains the mapping of ensembllink ids to GO terms
gene2GO <- list()

gene2GO[["Mm"]] <- by(annot.frame, annot.frame$ensembl_id,
                      function(x) as.character(x$go_id))
## define the universe of tested genes
gene2GO[["Mm"]] <- gene2GO[["Mm"]][names(gene2GO[["Mm"]])%in%
                                       rownames(Mm.1st.pass.model[[1]])]


go_Ef <- read.csv("data/12864_2014_6777_GO-annot.csv")

gene2GO[["Ef"]] <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))
gene2GO[["Ef"]] <- gene2GO[["Ef"]][names(gene2GO[["Ef"]])%in%
                                       rownames(Ef.1st.pass.model[[1]])]


#####  objects created here for further use:
## 1. annot.frame - a data frame containing coprehensive annotation information
## 2. gene2GO - a list linking cufflinks ids to GO terms for both mouse and Eimeria
