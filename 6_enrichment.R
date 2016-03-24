## still open question: How to implement the FCS(?) method
library(topGO)
library(xtable)

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
    cbind(cluster=2,type="oocysts", species="Ef"),
    cbind(cluster=5,type="sporozoites", species="Ef"),
    cbind(cluster=1,type="day7", species="Ef"),
    ## Mouse
    cbind(cluster=2,type="7dpiUp", species="Mm"),
    cbind(cluster=3,type="7dpiDown", species="Mm"),
    cbind(cluster=1,type="EarlyUp", species="Mm"),
    cbind(cluster=4,type="allInfDown", species="Mm")
    ))

apply(to.test, 1, function (x){
          clus = x[[1]]
          type = x[[2]]
          species = x[[3]]
          g2G = gene2GO[[species]]
          hcl = hcluster[[species]]
          file.path=("~/Ef_RNAseq/additional_files/tex/cluster")
          sapply(c("MF", "BP"), function (onto){
                     res <- TOGO.all.onto(onto, exp.universe[[species]],
                                          names(hcl[hcl%in%clus]),g2G)
                     file.detail= paste(clus, "_", onto, "_",
                         species, "_", type, ".tex", sep="")
                     capture.output(print(xtable(gene.table.topGO(res)),
                                          file=paste(file.path, file.detail, sep="")))
                 })
      })

