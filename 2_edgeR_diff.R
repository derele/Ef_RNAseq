#### DIFFERENTIAL GENE EXPRESSION ANALYSIS
#### OF RNA-seq DATA
## LIBRARIES
library(edgeR)
library(ggplot2)
library(GGally)

if(!exists("mouse.bg")){
    source("1_ballgown_import_diff.R")
}

raw.counts.4.bg <- function(bg){
exon.raw.count <- eexpr(bg, "ucount")
## the linkage data
e2t <- bg@indexes$e2t
t2g <- bg@indexes$t2g 
e2t2g <- merge(e2t, t2g)
all.count <- merge(e2t2g, exon.raw.count,
                   by.x = "e_id", by.y = 0)
## sum up for transcripts
transcript.raw.count <-
    do.call("rbind",
            by(all.count, all.count$t_id,
               function(x) colSums(x[, 4:ncol(x)])))
## sum up for genes
gene.raw.count <-
    do.call("rbind",
            by(all.count, all.count$g_id,
               function(x) colSums(x[, 4:ncol(x)])))
colnames(transcript.raw.count) <-
    as.character(pData(bg)$samples)
colnames(gene.raw.count) <-
    as.character(pData(bg)$samples)
return(list(exon.raw.count,
            transcript.raw.count,
            gene.raw.count))
}

mouse.RC <- raw.counts.4.bg(mouse.bg)

pdf("figures/rep_pairs_mouseG.pdf", width = 70, height=70)
ggpairs(mouse.RC[[3]])
dev.off()

## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
design <- model.matrix(~0+pData(mouse.bg)$grouped)

colnames(design)  <- gsub("pData\\(mouse.bg\\)\\$grouped", "",
                             colnames(design))
## just to have the standard point delimitor
colnames(design)  <- gsub("_", ".",
                             colnames(design))

keep <- rowSums(mouse.RC[[3]])>5000 ## filter: with 2000 as cutoff, bimodal distr. almost not visible
## tried 100, 1000, 5000. Latter is slightly better, but not much by visual inspection.
GM <- DGEList(mouse.RC[[3]][keep,], group = pData(mouse.bg)$grouped)

GM <- calcNormFactors(GM)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM <- estimateDisp(GM, design = design, robust = TRUE)

pdf("figures/rep_pairs_NORM_mouseG.pdf", width = 70, height=70)
ggpairs(cpm(GM))
dev.off()

pdf("figures/Mouse_genes_mds.pdf")
plotMDS(GM, )
dev.off()

pdf("figures/Mouse_genes_hclust.pdf")
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


