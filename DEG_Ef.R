## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

#setwd("/home/totta/Ef_RNAseq")

if(!exists("raw.counts.4.bg")){
	source("2_edgeR_diff.R")
        }

Ef.RC <- raw.counts.4.bg(Ef.bg)

pdf("figures/rep_pairs_EfG.pdf", width = 70, height=70)
ggpairs(Ef.RC[[3]])
dev.off()

## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
## For EIMERIA:
## need to rm day 0 from pData and design
design <- model.matrix(~0+pData(Ef.bg)$grouped)

colnames(design)  <- gsub("pData\\(Ef.bg\\)\\$grouped", "",
                             colnames(design))
## just to have the standard point delimitor
colnames(design)  <- gsub("_", ".",
                             colnames(design))

keep <- rowSums(Ef.RC[[3]])>100 ## filter lowly expressed
GM <- DGEList(Ef.RC[[3]][keep,], group = pData(Ef.bg)$grouped)

GM <- calcNormFactors(GM)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM <- estimateDisp(GM, design = design, robust = TRUE)

pdf("figures/rep_pairs_NORM_EfG.pdf", width = 70, height=70)
ggpairs(cpm(GM))
#dev.off()

pdf("figures/Ef_genes_mds.pdf")
plotMDS(GM, )
dev.off()

pdf("figures/Ef_genes_hclust.pdf")
plot(hclust(dist(t(cpm(GM)))))
dev.off()

fitM <- glmFit(GM, design) 

my.contrasts <- makeContrasts(
    ## just infections in NMRI
    N3vs0 = NMRI.1stInf.3dpi-NMRI.1stInf.0dpi,
    N5vs0 = NMRI.1stInf.5dpi-NMRI.1stInf.0dpi,
    N7vs0 = NMRI.1stInf.7dpi-NMRI.1stInf.0dpi,
    ## just infections in Black and Rag
    B5vs0 = C57BL6.1stInf.5dpi- C57BL6.0dpi ,
    R5vs0 = Rag.1stInf.5dpi - Rag.0dpi,
    ## differences in 1st vs 2nd infection
    N0.1stvsN0.2nd = NMRI.1stInf.0dpi - NMRI.2ndInf.0dpi,
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR2.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,
    ## differences in mouse strains
        NvsB=NMRI.1stInf.0dpi-C57BL6.0dpi, 
        NvsR=NMRI.1stInf.0dpi-Rag.0dpi,
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
    levels=design)

fitLRT <- glmLRT(fitM, contrast = my.contrasts)

ALL.top <- topTags(fitLRT, 1000000)

fitLRT.list <- lapply(colnames(my.contrasts),
                          function (x) glmLRT(fitM, contrast = my.contrasts[,x]))
names(fitLRT.list) <- colnames(my.contrasts)

top.list <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
names(top.list) <- colnames(my.contrasts)

gene.list <- lapply(top.list, function(x) {
                            rownames(x[x$FDR<0.01,])
                        })
