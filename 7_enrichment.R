library(topGO)
library(xtable)
library(AnnotationDbi)

if(!exists("hcluster")){
    source("5_clustering.R")
}
annotation.frame <- read.table("output_data/annotation_data.csv", sep=",")

## just to get the universe of genes which were tested at all:
Mm.DE.test <- read.table("output_data/Mm_DEtest.csv", sep=",")
Ef.DE.test <- read.table("output_data/Ef_DEtest.csv", sep=",")

gene2GO <- list()
gene2GO[["Mm"]] <- by(annot.frame, annot.frame$ensembl_id,
                      function(x) as.character(x$go_id))
## the universe of genes which were tested at all
gene2GO[["Mm"]] <- gene2GO[["Mm"]][names(gene2GO[["Mm"]])%in%
                                   unique(Mm.DE.test$gene)]


go_Ef <- read.csv("data/12864_2014_6777_GO-annot.csv")

gene2GO[["Ef"]] <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))
## the universe of genes which were tested at all
gene2GO[["Ef"]] <- gene2GO[["Ef"]][names(gene2GO[["Ef"]])%in%
                                   unique(Ef.DE.test$gene)]

exp.universe <- list()
exp.universe[["Mm"]] <- names(gene2GO[["Mm"]])
exp.universe[["Ef"]] <- names(gene2GO[["Ef"]])


hcluster <- list()
hcluster[["Mm"]] <- read.table("output_data/Mm_hclustered_cycle.csv", sep=",",
                               header=TRUE)
hcluster[["Ef"]] <- read.table("output_data/Ef_hclustered_cycle.csv", sep=",",
                               header=TRUE)

#### We search for enrichment of genes from clustering 
######## Clustering #############
to.test <- as.data.frame(rbind(
    ## Eimeria
    ## cbind(cluster=1,type="Day7Up_1", species="Ef"),
    ## cbind(cluster=2,type="Day7Up_2", species="Ef"),
    ## cbind(cluster=3,type="OocSporDown", species="Ef"),
    ## cbind(cluster=4,type="OocystsUp", species="Ef"),
    ## cbind(cluster=5,type="OoSpoDown_day5Up", species="Ef"),
    ## cbind(cluster=6,type="SporozoitesUp", species="Ef"),
    ## cbind(cluster=7,type="Day7OocyDown_SporoUP", species="Ef"),
    ## Mouse
    cbind(cluster=5,type="7dpiUp", species="Mm"),
    cbind(cluster=3,type="allInfDown", species="Mm"),
    cbind(cluster=4,type="CompInfDown", species="Mm")
    ))

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
    all$p.value <- gsub(" ?< ?", "", all$p.value)
    all$p.value <- as.numeric(all$p.value)
    all$adj.p <- round(p.adjust(all$p.value, method="BH"), 4)
    return(all[all$p.value<pval,])
}

## Subset BPMF.ll by [[1]][1[]][1[]] to get GO-terms (third bracket)
## per cluster (first No.) and either BP (1) or MF (2) (second
## number).
BPMF.ll <- apply(to.test, 1, function (x){
          clus = x[[1]]
          type = x[[2]]
          species = x[[3]]
          g2G = gene2GO[[species]] # gene2GO object created in 3_annotations script
          hcl = as.data.frame(hcluster[[species]])
          ## creates first part of filename for each cluster
          file.path=("~/Ef_RNAseq/Supplement/tex/cluster")
          BPMF.l <- lapply(c("MF", "BP"), function (onto){
              res <- TOGO.all.onto(onto, exp.universe[[species]],
                                   rownames(hcl)[hcl$Cluster%in%clus],g2G)
              file.detail= paste(clus, "_", onto, "_",
                                 species, "_", type, ".tex", sep="")
              capture.output(print(xtable(gene.table.topGO(res)),
                                   include.rownames=FALSE,
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





