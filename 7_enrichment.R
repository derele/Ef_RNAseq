library(topGO)
library(xtable)

annotation.frame <- read.table("output_data/annotation_data.csv", sep=",")

##  get the universe of genes which were tested at all:
Mm.DE.test <- read.table("output_data/Mm_DEtest.csv", sep=",")
Ef.DE.test <- read.table("output_data/Ef_DEtest.csv", sep=",")

gene2GO <- list()
gene2GO[["Mm"]] <- by(annotation.frame, annotation.frame$ensembl_id,
                      function(x) as.character(x$go_id))

## the universe of genes which were tested at all
gene2GO[["Mm"]] <- gene2GO[["Mm"]][names(gene2GO[["Mm"]])%in%
                                   unique(Mm.DE.test$gene)]

go_Ef <- read.csv("data/12864_2014_6777_GO-annot.csv")

gene2GO[["Ef"]] <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))
## the universe of genes which were tested at all
gene2GO[["Ef"]] <- gene2GO[["Ef"]][names(gene2GO[["Ef"]])%in%
                                   unique(Ef.DE.test$gene)]

##  get the universe of genes which were tested at all:
exp.universe <- list()
exp.universe[["Mm"]] <- names(gene2GO[["Mm"]])
exp.universe[["Ef"]] <- names(gene2GO[["Ef"]])

hcluster <- list()
hcluster[["Mm"]] <- read.table("output_data/Mm_hclustered_cycle.csv", sep=",",
                               header=TRUE)
hcluster[["Ef"]] <- read.table("output_data/Ef_hclustered_cycle.csv", sep=",",
                               header=TRUE)

set.from.DE <- function(test.df, contrast, FDR=0.01){
    target.df <- test.df[test.df$contrast%in%contrast, ]
    as.character(target.df[target.df$FDR<FDR, "gene"])
}


set.from.cluster <- function(hcluster, number){
    as.character(rownames(hcluster)[hcluster$Cluster%in%number])
}

#### We search for enrichment of genes from clustering 
######## Clustering #############
to.test <- list(
    ## Eimeria
    ## fill me in!!!
    ## Mouse
    list(set=set.from.cluster(hcluster[["Mm"]], 5),
         type="7dpiUp", species="Mm"),
    list(set=set.from.cluster(hcluster[["Mm"]], 3),
          type="allInfDown", species="Mm"),
    list (set=set.from.cluster(hcluster[["Mm"]], 4),
          type="CompInfDown", species="Mm"),
    list (set=set.from.DE(Mm.DE.test, "N5vsN7"),
          type="DE_N5vsN7", species="Mm")
    )

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
BPMF.ll <- lapply(to.test, function (x){
          set = x[[1]]
          type = x[[2]]
          species = x[[3]]
          g2G = gene2GO[[species]] # gene2GO object created in 3_annotations script
          ## creates first part of filename for each cluster
          file.path=("Supplement/tex/cluster")
          BPMF.l <- lapply(c("MF", "BP"), function (onto){
              res <- TOGO.all.onto(onto, exp.universe[[species]],
                                   set, g2G)
              file.detail= paste(onto, "_",
                                 species, "_", type, ".tex", sep="")
              capture.output(print(xtable(gene.table.topGO(res)),
                                   include.rownames=FALSE,
                                   file=paste(file.path, file.detail, sep="")))
              return(gene.table.topGO(res))
          })
          names(BPMF.l) <- c("MF", "BP")
          return(type=BPMF.l)
})

names(BPMF.ll) <- unlist(lapply(to.test, "[[", 2))

### NOT tested from here on...
## Short how-to:

## 1. Add the gene set as first element in the list to.test

##2. Add a short (one word) description of the gene set as second
##element

## 3. Run to get genes contributing to BP and MF enrichments in
## cluster/DE Genes of interest

## Get all GO-terms from one cluster
clusterMF <- BPMF.ll[[4]][[2]][[1]]
clusterBP <- BPMF.ll[[4]][[1]][[1]]
clusterBPMF <- append(clusterBP, clusterMF) # appends unique GO:terms

## Get ancestor annotations
frame.data <-  data.frame(cbind(GOIDs=as.character(go_Ef$go),
                                evi="ND", genes=as.character(go_Ef$gene)))

frame=GOFrame(frame.data, organism="Eimeria falciformis")

allFrame=GOAllFrame(frame)

go_Ef_all <- getGOFrameData(allFrame)

## Get genes names contributing to specific GO-term in specific cluster
genes.from.subset <-
    as.character(unique(go_Ef_all[go_Ef_all$go_id%in%clusterBPMF, 3])) 
sign.genes <-
    genes.from.subset[genes.from.subset%in%names(hcluster[["Ef"]])
                      [hcluster[["Ef"]]%in%4]] 

## ALL  genes in each cluster (adj. cluster number)
all.GOgenes.in.cluster <-
    names(hcluster[["Ef"]])[hcluster[["Ef"]]%in%1]





