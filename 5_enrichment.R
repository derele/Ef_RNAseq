library(topGO)

if(!exists("chalGids")){
    source("1_ballgown_import_diff.R")
}

if(!exists("gene.list")){
    source("2_edgeR_diff.R")
}

if(!exists("gene2GO")){
    source("3_annotations.R")
}


chalGids <- get.annotation.for.xloc(chal_genes)[[2]]

MF.chal <- TOGO.all.onto("MF", names(gene2GO), 
                         chalGids, gene2GO)
gene.table.topGO(MF.chal)

BP.chal <- TOGO.all.onto("BP", names(gene2GO), 
                         chalGids , gene2GO)
head(gene.table.topGO(BP.chal), n=20)



chalGidsI <- get.annotation.for.xloc(chal_genesI)[[2]]

MF.chalI <- TOGO.all.onto("MF", names(gene2GO), 
                         chalGidsI, gene2GO)
gene.table.topGO(MF.chalI)

BP.chalI <- TOGO.all.onto("BP", names(gene2GO), 
                          chalGidsI , gene2GO)
head(gene.table.topGO(BP.chalI, 1), n=20)



strainGids5 <- get.annotation.for.xloc(strain_genes5)[[2]]

MF.strain5 <- TOGO.all.onto("MF", names(gene2GO), 
                            strainGids5, gene2GO)
gene.table.topGO(MF.strain5)

BP.strain5 <- TOGO.all.onto("BP", names(gene2GO), 
                            strainGids5, gene2GO)
head(gene.table.topGO(BP.strain5), n=20)


strainGidsI <- get.annotation.for.xloc(strain_genesI)[[2]]

MF.strainI <- TOGO.all.onto("MF", names(gene2GO), 
                            strainGidsI, gene2GO)
gene.table.topGO(MF.strainI)

BP.strainI <- TOGO.all.onto("BP", names(gene2GO), 
                            strainGidsI, gene2GO)
head(gene.table.topGO(BP.strainI), n=20)



