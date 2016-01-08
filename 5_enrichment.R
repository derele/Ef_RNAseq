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

## the universe of genes which were tested at all:
exp.universe <- get.annotation.for.xloc(gene.list[[2]])[[2]]


## The CHALLENGE in NMRI DAY 7 
N7.1stvsN7.2nd.GID <-
    get.annotation.for.xloc(gene.list[["N7.1stvsN7.2nd"]])[[2]]

MF.N7.chal <- TOGO.all.onto("MF", exp.universe,
                            N7.1stvsN7.2nd.GID, gene2GO)
gene.table.topGO(MF.N7.chal)

BP.N7.chal <- TOGO.all.onto("BP", exp.universe, 
                            N7.1stvsN7.2nd.GID, gene2GO)
gene.table.topGO(BP.N7.chal)

topKEGG(kegga(N7.1stvsN7.2nd.GID, species = "Mm"), n=30)
topGO(goana(N7.1stvsN7.2nd.GID, species = "Mm"), n=30)


## the genes not detected as DE  on the array but in RNA seq:
Array.universe <-
    get.annotation.for.xloc(RNAseq.Array.logFC$Row.names)[[2]]

A.R.fail <- get.annotation.for.xloc(A.low.R.high)[[2]]

AR.fail.MF <- TOGO.all.onto("MF", Array.universe,
                            A.R.fail, gene2GO)
gene.table.topGO(AR.fail.MF)

AR.fail.BP <- TOGO.all.onto("BP", Array.universe,
                            A.R.fail, gene2GO)
gene.table.topGO(AR.fail.MF)

topKEGG(kegga(A.R.fail, species = "Mm"), n=30)
topGO(goana(A.R.fail, species = "Mm"), n=30)






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



