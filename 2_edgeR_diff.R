#### DIFFERENTIAL GENE EXPRESSION ANALYSIS

## needed from data import script are two raw count objects 

if(!exists("Mm.RC")|!exists("Ef.RC")){
    source("1_ballgown_import.R")
}

library(edgeR)
library(GGally)

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
Mm.design <- model.matrix(~0+pData(Mm.bg)$grouped)

colnames(Mm.design)  <- gsub("pData\\(Mm.bg\\)\\$grouped", "",
                             colnames(Mm.design))
## just to have the standard point delimitor
colnames(Mm.design)  <- gsub("_", ".",
                             colnames(Mm.design))

 ## filter: with 3000 as cutoff, bimodal distr. almost not visible
keep <- rowSums(Mm.RC[[3]])>3000

DGEList.Mm <- DGEList(Mm.RC[[3]][keep,], group = pData(Mm.bg)$grouped)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
DGEList.Mm <- estimateDisp(DGEList.Mm, design = Mm.design, robust = TRUE)

fit.Mm <- glmFit(DGEList.Mm, Mm.design) 

Mm.contrasts <- makeContrasts(
    ## just infections in NMRI
    N3vs0 = NMRI.1stInf.3dpi-NMRI.2ndInf.0dpi,
    N5vs0 = NMRI.1stInf.5dpi-NMRI.2ndInf.0dpi,
    N7vs0 = NMRI.1stInf.7dpi-NMRI.2ndInf.0dpi,
    ## just infections in Black and Rag
    B5vs0 = C57BL6.1stInf.5dpi- C57BL6.0dpi ,
    R5vs0 = Rag.1stInf.5dpi - Rag.0dpi,
    ## differences in 1st vs 2nd infection
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR2.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,
    ## differences in mouse strains
        RvsB=Rag.0dpi-C57BL6.0dpi,
    ## GO on here:
    ## differences in infection over time
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi,
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi,
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi,
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi),
    levels=Mm.design)

Mm.fitLRT <- glmLRT(fit.Mm, contrast = Mm.contrasts)

ALL.top.Mm <- topTags(Mm.fitLRT, 1000000)

Mm.fitLRT.list <- lapply(colnames(Mm.contrasts),
                         function (x) glmLRT(fit.Mm,
                                             contrast = Mm.contrasts[,x]))
names(Mm.fitLRT.list) <- colnames(Mm.contrasts)

Mm.top.list <- lapply(Mm.fitLRT.list, function (x){
                          topTags(x, 1000000)[[1]]
                       })
names(Mm.top.list) <- colnames(Mm.contrasts)

Mm.gene.list <- lapply(Mm.top.list, function(x) {
                           rownames(x[x$FDR<0.01,])
                       })

