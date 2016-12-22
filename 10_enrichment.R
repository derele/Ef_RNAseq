library(topGO)
library(xtable)
library(gridExtra) #  for color theme function
library(grid)      # for color theme function
library(reshape)

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

interpro_Ef <- read.csv("data/12864_2014_6777_interpro_EfaB.csv")


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

test.GO.control.list <- function (x){
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
}

BPMF.ll <- lapply(to.test, test.GO.control.list)

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

Ef.tested.universe <- unique(Ef.DE.test$gene)
Mm.tested.universe <- unique(Mm.DE.test$gene)

IPR <- read.delim("data/ipr_proteins.fa.tsv", header=FALSE,
                  as.is=TRUE)

SigP_euk <- IPR$V1[IPR$V4%in%c("SignalP_EUK")]

SigP <- IPR$V1[IPR$V4%in%c("SignalP_EUK",
                           "SignalP_GRAM_NEGATIVE",
                           "SignalP_GRAM_POSITIVE")]

SigTMHMM <- IPR$V1[IPR$V4%in%c("TMHMM")]

get.unique.SigP.genes <- function(x){
    u <- unname(unlist(sapply(x, strsplit, "\\|")))
    u <- unique(gsub(".t\\d+", "", u))
    gsub("NODE_", "EfaB_", u)
}

SigP_euk <- get.unique.SigP.genes(SigP_euk)
SigP <- get.unique.SigP.genes(SigP)
SigTMHMM <- get.unique.SigP.genes(SigTMHMM)

cluster.p.SigP <- lapply(unique(hcluster[["Ef"]]$Cluster), function(x){
    ft <- fisher.test(Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], x),
                      Ef.tested.universe %in% SigP)
    list(ft$estimate, ft$p.value)
})


cluster.p.SigP_euk <- lapply(unique(hcluster[["Ef"]]$Cluster), function(x){
    ft <- fisher.test(Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], x),
                   Ef.tested.universe %in% SigP_euk)
     list(ft$estimate, ft$p.value)
})


cluster.p.TMHMM <- lapply(unique(hcluster[["Ef"]]$Cluster), function(x){
    ft <- fisher.test(Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], x),
                      Ef.tested.universe %in% SigTMHMM)
    list(ft$estimate, ft$p.value)
})

### Cluster 5 (oocysts) enriched for transmembrane signal
table(Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], 5),
      Ef.tested.universe %in% SigTMHMM)

cluster.5.Membrane.genes <- Ef.tested.universe[
    Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], 5) &
    Ef.tested.universe %in% SigTMHMM]

write.csv(cluster.5.Membrane.genes, "data/Cluster_5_Membrane.txt")

SigClus <- data.frame(rbind(do.call(rbind, cluster.p.SigP),
                            do.call(rbind, cluster.p.SigP_euk),
                            do.call(rbind, cluster.p.TMHMM)))

names(SigClus) <- c("odds.ratio", "p.value")

SigClus$cluster <- 1:7

SigClus$test <- rep(c("SigP", "SigP_euk", "TMHMM"), each=7)

SigClus$adj.p <- p.adjust(SigClus$p.value, method="BH")



############### R&B
## load the full products in RnB.final
load("/SAN/Eimeria_Totta/RnB_Prod_1478263755.Rdata")

is.zero <- (RnB.final==0)

is.zero.all.cols <- apply(is.zero, 2, any)
Ef.genes.interacting <- colnames(RnB.final)[is.zero.all.cols] 

is.zero.all.rows <- apply(is.zero, 1, any)
Mm.genes.interacting <- rownames(RnB.final)[is.zero.all.rows]

## interaction clusters ... R&B 
Ef.interA.p.Cluster <- lapply(unique(hcluster[["Ef"]]$Cluster), function(x){
    cluster.set <- set.from.cluster(hcluster[["Ef"]], x)
    n.cluster <- length(cluster.set)
    N.genes.interacting <- length(cluster.set[cluster.set %in% Ef.genes.interacting])
    perc.genes.interacting <- (N.genes.interacting/n.cluster)*100
    ft <- fisher.test(Ef.tested.universe %in% set.from.cluster(hcluster[["Ef"]], x),
                      Ef.tested.universe %in% Ef.genes.interacting)
    list(length(cluster.set), perc.genes.interacting, ft$estimate, ft$p.value)
})

Ef.interA.p.Cluster <- data.frame(cbind(do.call(rbind, Ef.interA.p.Cluster),
                                        Cluster=paste("E.falciformis", 1:7)))
                                  
names(Ef.interA.p.Cluster) <- c("N cluster", "% interacting",
                                "odds ratio", "p.value", "Cluster")
Ef.interA.p.Cluster$FDR <- p.adjust(Ef.interA.p.Cluster$p.value, method="BH")

## this table could be reported, but better just the tow significant values
Ef.interA.p.Cluster

Mm.interA.p.Cluster <- lapply(unique(hcluster[["Mm"]]$Cluster), function(x){
    cluster.set <- set.from.cluster(hcluster[["Mm"]], x)
    n.cluster <- length(cluster.set)
    N.genes.interacting <- length(cluster.set[cluster.set %in% Mm.genes.interacting])
    perc.genes.interacting <- (N.genes.interacting/n.cluster)*100
    ft <- fisher.test(Mm.tested.universe %in% set.from.cluster(hcluster[["Mm"]], x),
                      Mm.tested.universe %in% Mm.genes.interacting)
    list(length(cluster.set), perc.genes.interacting, ft$estimate, ft$p.value)
})

Mm.interA.p.Cluster <- data.frame(cbind(do.call(rbind, Mm.interA.p.Cluster),
                                        Cluster=paste("Mouse", 1:7)))
                                  
names(Mm.interA.p.Cluster) <- c("N cluster", "% interacting",
                                "odds ratio", "p.value", "Cluster")
Mm.interA.p.Cluster$FDR <- p.adjust(Mm.interA.p.Cluster$p.value, method="BH")

## this table could be reported, but better just the tow significant values
Mm.interA.p.Cluster

Clust.Enrich <- rbind(Mm.interA.p.Cluster, Ef.interA.p.Cluster)

table.clust.tex <- xtable(Clust.Enrich, digits=c(NA, 0, 2, 2, -2, NA, -2))

print(table.clust.tex,
      type = "html", file = "tables/Table4_ISIGM_Cluster.html", include.rownames = F,
      format.args = list(big.mark = ",", decimal.mark = "."))

## here the quantitative insight into this:
RnB.cluster.scores <-
    lapply(unique(hcluster[["Mm"]]$Cluster), function(x){
        lapply(unique(hcluster[["Ef"]]$Cluster), function(y){
            Clus.set.Ef <- set.from.cluster(hcluster[["Ef"]], y)
            Clus.set.Mm <- set.from.cluster(hcluster[["Mm"]], x)
            RnB.clus <- RnB.final[Clus.set.Mm, Clus.set.Ef]
            as.vector(RnB.clus)
        })
    })

isigem.cluster.in.cluster <- melt(RnB.cluster.scores)
names(isigem.cluster.in.cluster) <- c("ISIGEM.Score", "EfExpCluster", "MmExpCluster")
isigem.cluster.in.cluster$EfExpCluster <-
    as.factor(paste0("EfCluster",
                     isigem.cluster.in.cluster$EfExpCluster))
isigem.cluster.in.cluster$MmExpCluster <-
    as.factor(paste0("MmCluster",
                     isigem.cluster.in.cluster$MmExpCluster))

## plotting the quantitative view on ISIGM scores
devSVG("figures/Figure5d_IsigemClusterInCluster.svg")
ggplot(isigem.cluster.in.cluster, aes(ISIGEM.Score, ..count.., colour=MmExpCluster)) +
    geom_density() +
    facet_wrap(~EfExpCluster) +
    theme_bw()
dev.off()

fisher.test(Ef.tested.universe %in% SigTMHMM,
            Ef.tested.universe %in% Ef.genes.interacting)

fisher.test(Ef.tested.universe %in% SigP_euk,
            Ef.tested.universe %in% Ef.genes.interacting)

fisher.test(Ef.tested.universe %in% SigP,
            Ef.tested.universe %in% Ef.genes.interacting)
## all NS


Clus.7dpi <- rownames(hcluster[["Ef"]])[hcluster[["Ef"]]$Cluster==2]

Inter.DE <- Ef.genes.interacting[Ef.genes.interacting%in%Clus.7dpi]


## underrepresentation of genes in interacting
fisher.test(Ef.tested.universe %in% SigP,
            Ef.tested.universe %in% Inter.DE)

to.test.2 <- list(
    ## Eimeria
    list(set=Inter.DE,  
         type = "Interacting7dpi", species = "Ef"),
    list(set=Ef.genes.interacting,  
         type = "Interacting", species = "Ef")
)

BPMF.inter.ll <- lapply(to.test.2, test.GO.control.list)

names(BPMF.inter.ll) <- unlist(lapply(to.test.2, "[[", 2))
