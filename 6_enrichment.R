## still open question: How to implement the FCS(?) method
# Check checkpoint package for package control.
if(!require(topGO)) biocLite("topGO") # Imports package if user does not have it
if(!require(xtable)) biocLite("xtable") # Imports package if user does not have it
if(!require(AnnotationDbi)) biocLite("AnnotationDbi") # Imports package if user does not have it

library(topGO)
library(xtable)
library(AnnotationDbi)

if(!exists("hcluster")){
    source("5_clustering.R")
}

if(!exists("gene2GO")){
    source("3_annotations.R")
}

## the universe of genes which were tested at all:
exp.universe <- list()
exp.universe[["Mm"]] <- names(gene2GO[["Mm"]])
exp.universe[["Ef"]] <- names(gene2GO[["Ef"]])


#### We search for enrichment of genes from clustering 
######## Clustering #############

to.test <- as.data.frame(rbind(
    ## Eimeria
    cbind(cluster=1,type="Day7Up_1", species="Ef"),
    cbind(cluster=2,type="Day7Up_2", species="Ef"),
    cbind(cluster=3,type="OocSporDown", species="Ef"),
    cbind(cluster=4,type="OocystsUp", species="Ef"),
    cbind(cluster=5,type="OoSpoDown_day5Up", species="Ef"),
    cbind(cluster=6,type="SporozoitesUp", species="Ef"),
    cbind(cluster=7,type="Day7OocyDown_SporoUP", species="Ef"),
    ## Mouse
    cbind(cluster=2,type="7dpiUp", species="Mm"),
    cbind(cluster=3,type="7dpiDown", species="Mm"),
    cbind(cluster=1,type="EarlyUp", species="Mm"),
    cbind(cluster=4,type="allInfDown", species="Mm")
    ))


## Subset BPMF.ll by [[1]][1[]][1[]] to get GO-terms (third bracket)  per cluster (first No.) and 
## either BP (1) or MF (2) (second number).
BPMF.ll <- apply(to.test, 1, function (x){
          clus = x[[1]]
          type = x[[2]]
          species = x[[3]]
          g2G = gene2GO[[species]] # gene2GO object created in 3_annotations script
          hcl = hcluster[[species]] # hcluster object created in 5_clustering script
          ## creates first part of filename for each cluster
          file.path=("~/Ef_RNAseq/additional_files/tex/cluster")
          BPMF.l <- lapply(c("MF", "BP"), function (onto){
                       res <- TOGO.all.onto(onto, exp.universe[[species]],
                                            names(hcl[hcl%in%clus]),g2G)
                      file.detail= paste(clus, "_", onto, "_",
                           species, "_", type, ".tex", sep="")
                      capture.output(print(xtable(gene.table.topGO(res)), include.rownames=FALSE,
                                            file=paste(file.path, file.detail, sep="")))
                      return(gene.table.topGO(res))
                 })
          return(BPMF.l)
      })
 
## 1. Short how-to: change number in first bracket according to cluster of interest
## 2. Adjust cluster number where sign.genes is created
## 3. Run to get genes contributing to BP and MF enrichments in cluster of interest

## Get all GO-terms from one cluster
clusterMF <- BPMF.ll[[4]][[2]][[1]]
clusterBP <- BPMF.ll[[4]][[1]][[1]]
clusterBPMF <- append(clusterBP, clusterMF) # appends unique GO:terms

## Get ancestor annotations
frame.data <- data.frame(cbind(GOIDs=as.character(go_Ef$go), evi="ND", genes=as.character(go_Ef$gene)))
frame=GOFrame(frame.data, organism="Eimeria falciformis")
allFrame=GOAllFrame(frame)
go_Ef_all <- getGOFrameData(allFrame)

## Get genes names contributing to specific GO-term in specific cluster
genes.from.subset <- as.character(unique(go_Ef_all[go_Ef_all$go_id%in%clusterBPMF, 3])) 
sign.genes <- genes.from.subset[genes.from.subset%in%names(hcluster[["Ef"]])[hcluster[["Ef"]]%in%4]] 
sign.genes
## ALL  genes in each cluster (adj. cluster number)
all.GOgenes.in.cluster <- names(hcluster[["Ef"]])[hcluster[["Ef"]]%in%1]





