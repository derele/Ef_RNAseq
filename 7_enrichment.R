library(topGO)
library(xtable)
library(gridExtra) #  for color theme function
library(grid)      # for color theme function

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

## EIMERIA
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
    list(set=set.from.cluster(hcluster[["Ef"]], 1),  # no enriched MF terms
         type = "cluster1ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 5),
         type = "cluster5ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 7),
         type = "cluster7ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 2),
         type = "cluster2ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 4),
          type = "cluster4ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 6),
         type = "cluster6ef", species = "Ef"),
    list(set=set.from.cluster(hcluster[["Ef"]], 3),   # no enriched BP terms
         type = "cluster3ef", species = "Ef"),
    ## Mouse
    list(set=set.from.cluster(hcluster[["Mm"]], 6),
         type="cluster6mm", species="Mm"),
    list(set=set.from.cluster(hcluster[["Mm"]], 5),
         type="cluster5mm", species="Mm"),
    list(set=set.from.cluster(hcluster[["Mm"]], 3),
          type="cluster3mm", species="Mm"),
    list (set=set.from.cluster(hcluster[["Mm"]], 4),
          type="cluster4mm", species="Mm"),
    list (set=set.from.cluster(hcluster[["Mm"]], 1),
          type="cluster1mm", species="Mm"),
    list (set=set.from.cluster(hcluster[["Mm"]], 2),
          type="cluster2mm", species="Mm"),
    list (set=set.from.cluster(hcluster[["Mm"]], 7),
          type="cluster7mm", species="Mm"),
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

## Subset BPMF.ll by e.g. BPMF.ll$clusterxx1$BP (another $Term will give only annotated terms)
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

##colors for table
myt <- ttheme_default(
  base_size = 18,
  padding = unit(c(2, 6), "mm"),
  # Use hjust and x to left justify the text
  # Alternate the row fill colours
  core = list(fg_params=list(col="dark green"), #, hjust = 1, x=1),
              bg_params=list(fill=c("white", "light gray"))),
  
  # Change column header to white text and red background
  colhead = list(fg_params=list(col="dark green"),
                 bg_params=list(fill="gray"))
)

## Change filename and part of BPMF.ll object to export tables 
## (adjust size of PDF for better readability)
pdf("~/Ef_RNAseq/Supplement/SI_4_GO-enrichments-per-cluster/SI_GOMmcluster2_mf.pdf", width = 24, height = 10)
grid.table(data.frame(BPMF.ll$cluster2mm$BP),
                               theme = myt,
                               rows = NULL,
                               cols = c("GO id", "Term", 
                                        "Annot. genes", 
                                        "Sign. genes", 
                                        "Expected", "P-value", "Adj. P-value"))
dev.off()


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
get.Ef.ancestors <- function(){
frame.data <-  data.frame(cbind(GOIDs=as.character(go_Ef$go),
                                evi="ND", genes=as.character(go_Ef$gene)))
frame <- GOFrame(frame.data, organism="Eimeria falciformis")
allFrame <- GOAllFrame(frame)
getGOFrameData(allFrame)
}

go_Ef_all <- get.Ef.ancestors()


get.Mm.ancestors <- function(){
  frame.data <-  data.frame(cbind(GOIDs=as.character(annotation.frame$go_id),
                                  evi="ND", 
                                  genes=as.character(annotation.frame$ensembl_id)))
  frame <- GOFrame(frame.data, organism="Mus musculus")
  allFrame <- GOAllFrame(frame)
  getGOFrameData(allFrame)
}

go_Mm_all <- get.Mm.ancestors()


## helper functions
get.gene.4.go.Ef <- function (goterm){
  y <- go_Ef_all[go_Ef_all$go_id%in%goterm, "gene_id"]
  unique(y[y%in%exp.universe[["Ef"]]])
}

get.gene.4.go.Mm <- function (goterm){
  y <- go_Mm_all[go_Mm_all$go_id%in%goterm, "gene_id"]
  unique(y[y%in%exp.universe[["Mm"]]])
}

## use this for both mouse and eimeria; speciffy org if not Mm
get.genes.4.heatclus <- function(org="Mm", clust){
  rownames(hcluster[[org]])[hcluster[[org]]$Cluster==clust]
}

## execute function to get genes responsible for specific GO-term enrichments
get.go.set.Ef <- function(GO, clus){
  all <- get.gene.4.go.Ef(GO)
  all[all%in%get.genes.4.heatclus("Ef", clus)]
}

get.go.set.Mm <- function(GO, clus){
  all <- get.gene.4.go.Mm(GO)
  all[all%in%get.genes.4.heatclus("Mm", clus)]
}

