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

pdf("figures/rep_pairs_mouseG.pdf", width = 100, height=100)
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

GM <- DGEList(mouse.RC[[3]], group = pData(mouse.bg)$grouped)
GM <- calcNormFactors(GM)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM <- estimateDisp(GM, design = design, robust = TRUE)

pdf("figures/Mouse_genes_mds.pdf")
plotMDS(GM)
dev.off()

pdf("figures/Mouse_genes_hclust.svg")
plot(hclust(dist(t(cpm(GM)))))
dev.off()

fitM <- glmFit(GM, design) 
fitLRT.list <- lapply(1:6, function (i) glmLRT(fitM, coef=i))

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
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    ## differences in mouse strains
    NvsB=NMRI.1stInf.0dpi-C57BL6.0dpi, 
    NvsR=NMRI.1stInf.0dpi-Rag.0dpi,
    RvsB=Rag.0dpi-C57BL6.0dpi,
    ## GO on here:
    ## differences in infection over time
    ## differences in 1st vs. 2nd depending on mouse strain
    levels=design)

fitLRT <- glmLRT(fitM, contrast = my.contrasts)

ALL.top <- topTags(fitLRT, 1000000)[[1]]

fitLRT.list <- lapply(colnames(my.contrasts),
                          function (x) glmLRT(fitM, contrast = my.contrasts[,x]))
names(fitLRT.list) <- colnames(my.contrasts)

top.list <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
names(top.list) <- colnames(my.contrasts)

gene.list <- lapply(top.list, function(x) {
                            rownames(x[x$FDR<0.05,])
                        })


