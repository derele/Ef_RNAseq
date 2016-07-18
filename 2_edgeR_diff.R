#### DIFFERENTIAL GENE EXPRESSION ANALYSIS

## needed from data import script are two raw count objects 

#source("https://bioconductor.org/biocLite.R")
if(!exists("Mm.RC")|!exists("Ef.RC")){
    source("1_ballgown_import.R")
}


library(statmod)
library(edgeR)
library(GGally)
library(RUVSeq)
library(reshape)

## Implement the sample exclusion (see A_data_curation.R)
Mm.bg <- subset(Mm.bg,
                !pData(Mm.bg)$samples%in%c("NMRI_2ndInf_3dpi_rep1",
                                           "NMRI_2ndInf_5dpi_rep2",
                                           "NMRI_1stInf_0dpi_rep1"),
                genomesubset=FALSE)
Mm.RC <- raw.counts.4.bg(Mm.bg)

Ef.bg <- subset(Ef.bg,
                !pData(Ef.bg)$samples%in%c("NMRI_2ndInf_3dpi_rep1",
                                           "NMRI_2ndInf_5dpi_rep2",
                                           "NMRI_1stInf_0dpi_rep1"),
                genomesubset=FALSE)
Ef.RC <- raw.counts.4.bg(Ef.bg)

## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
Ef.design <- model.matrix(~0+pData(Ef.bg)$grouped)
colnames(Ef.design)  <- gsub("pData\\(Ef.bg\\)\\$grouped", "",
                             colnames(Ef.design))
## just to have the standard point delimitor
colnames(Ef.design)  <- gsub("_", ".",
                             colnames(Ef.design))

Ef.contrasts <- makeContrasts(
    ## just different stages in NMRI
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi, #1
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi,
    N3vsSpo = NMRI.1stInf.3dpi-NMRI.sporozoites,
    N3vsOoc = NMRI.1stInf.3dpi-NMRI.oocysts,
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi,
    N5vsSp = NMRI.1stInf.5dpi-NMRI.sporozoites,
    N5vsOo = NMRI.1stInf.5dpi-NMRI.oocysts,
    N7vsSp = NMRI.1stInf.7dpi-NMRI.sporozoites,
    N7vsOo = NMRI.1stInf.7dpi-NMRI.oocysts,
    SpvsOo=  NMRI.sporozoites-NMRI.oocysts,	#10
    ## differences in 1st vs 2nd infection
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,	#11
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR5.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,	#15
    ## differences in mouse strains
    NvsB=NMRI.1stInf.5dpi-C57BL6.1stInf.5dpi,	#16
    NvsR=NMRI.1stInf.5dpi-Rag.1stInf.5dpi,
    BvsR=C57BL6.1stInf.5dpi-Rag.1stInf.5dpi,	#18
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),	#19
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi),	#21
    levels=Ef.design)


## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
Mm.design <- model.matrix(~0+pData(Mm.bg)$grouped)

colnames(Mm.design)  <- gsub("pData\\(Mm.bg\\)\\$grouped", "",
                             colnames(Mm.design))
## just to have the standard point delimitor
colnames(Mm.design)  <- gsub("_", ".",
                             colnames(Mm.design))

Mm.contrasts <- makeContrasts(
    ## just infections in NMRI !!! As we not have uninfected control
    ## using 2nd infected control here
    N3vs0 = NMRI.1stInf.3dpi-NMRI.2ndInf.0dpi, #1
    N5vs0 = NMRI.1stInf.5dpi-NMRI.2ndInf.0dpi, #2
    N7vs0 = NMRI.1stInf.7dpi-NMRI.2ndInf.0dpi, #3
    ## just infections in Black and Rag
    B5vs0 = C57BL6.1stInf.5dpi- C57BL6.0dpi , #4
    R5vs0 = Rag.1stInf.5dpi - Rag.0dpi,		#5
    ## differences in 1st vs 2nd infection
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi, #6
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR5.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,	#10
    ## differences in mouse strains
        RvsB=Rag.0dpi-C57BL6.0dpi,
	#NvsB=NMRI.0dpi-C57BL6.0dpi,
    ## differences in infection over time
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi,
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi,
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi,		#14
    ## what about differences in the reaction to first infection in
    ## different mouse strains??
    ##
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),		#15
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi), 		#16
    levels=Mm.design)


get.my.models <- function (RC, cutoff, group,
                           design, contrasts,
                           norm="upperquartile"){
    keep <- rowSums(RC)>cutoff
    DGEList <- DGEList(RC[keep,])
    if(!is.null(norm)){
        DGEList <- calcNormFactors(DGEList, method=norm )
    }
    ## robust estimator https://support.bioconductor.org/p/65558/
    DGEList <- estimateDisp(DGEList, design = design, robust = TRUE)
    fit <- glmFit(DGEList, design) 
    fitLRT <- glmLRT(fit, contrast = contrasts)
    ALL.top <- topTags(fitLRT, 1000000)
    fitLRT.list <- lapply(colnames(contrasts),
                          function (x) glmLRT(fit,
                                              contrast = contrasts[,x]))
    names(fitLRT.list) <- colnames(contrasts)
    top.list <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
    names(top.list) <- colnames(contrasts)
    gene.list <- lapply(top.list, function(x) {
                            rownames(x[x$FDR<0.01,]) ######### CHANGE FDR!!!!!!!!!!
                        })
    return(list(ALL.top, top.list, gene.list, DGEList))
}

## filter: with 3000 as cutoff, bimodal distr. almost not visible (see
## A_data_curation.R)
Mm.1st.pass.model <- get.my.models(Mm.RC[[3]], cutoff=3000,
                                   group = pData(Mm.bg)$grouped,
                                   contrasts=Mm.contrasts,
                                   design=Mm.design,
                                   norm="upperquartile")

Ef.1st.pass.model <- get.my.models(Ef.RC[[3]], cutoff=100,
                                   group = pData(Ef.bg)$grouped,
                                   contrasts=Ef.contrasts,
                                   design=Ef.design,
                                   norm="upperquartile")


get.RUVed.data <- function(RC, cutoff, group){
    ## use the expression set and the filtering decided on before
    keep <- rowSums(RC)>cutoff
    RUVset <- newSeqExpressionSet(RC[keep,])
    RUVset <- betweenLaneNormalization(RUVset, which="upper")
    RUVs(RUVset, rownames(RUVset), k=1, makeGroups(group))
}

Mm.RUVset.groups <- get.RUVed.data(Mm.RC[[3]], cutoff = 3000,
                                   pData(Mm.bg)$grouped)

Ef.RUVset.groups <- get.RUVed.data(Ef.RC[[3]], cutoff = 100,
                                   pData(Ef.bg)$grouped)

Mm.RUVg.model <- get.my.models(normCounts(Mm.RUVset.groups),
                               cutoff=3000, ## no further cut off
                               group = pData(Mm.bg)$grouped,
                               contrasts=Mm.contrasts,
                               design=Mm.design, norm=NULL)

Ef.RUVg.model <- get.my.models(normCounts(Ef.RUVset.groups),
                               cutoff=100, ## no further cut off
                               group = pData(Ef.bg)$grouped,
                               contrasts=Ef.contrasts,
                               design=Ef.design, norm=NULL)

## how many genes are regulated:
cbind(melt(lapply(Mm.1st.pass.model[[3]], length)),
      melt(lapply(Mm.RUVg.model[[3]], length)))[,c(1,3,4)]

cbind(melt(lapply(Ef.1st.pass.model[[3]], length)),
      melt(lapply(Ef.RUVg.model[[3]], length)), by = "L1")[,c(1,3,4)]
      

