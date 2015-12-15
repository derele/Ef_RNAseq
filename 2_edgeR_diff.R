#### DIFFERENTIAL GENE EXPRESSION ANALYSIS
#### OF RNA-seq DATA
## LIBRARIES
library(edgeR)
library(ggplot2)


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


pdf("figures/rep_pairs_mouseG.pdf", width = 30, height=30)
ggpairs(mouse.RC[[3]])
dev.off()

design <- model.matrix(~pData(mouse.bg)$infection +
                           pData(mouse.bg)$mouse.strain +
                               pData(mouse.bg)$timepoint)

colnames(design)  <- gsub("pData\\(mouse.bg\\)\\$", "",
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

list.of.top.genesM <- lapply(fitLRT.list, topTags, 1000000)

names(list.of.top.genesM) <- unlist(lapply(list.of.top.genesM, 
                                    function(x) x@.Data[[3]]))

list.of.diff.genes <- lapply(list.of.top.genesM, function(x) {
                                 d  <- x[[1]]
                                 rownames(d[d$FDR<0.05,])
                             })

## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
designALT <- model.matrix(~0+pData(mouse.bg)$grouped)

colnames(designALT)  <- gsub("pData\\(mouse.bg\\)\\$grouped", "",
                             colnames(designALT))
## just to have the standard point delimitor
colnames(designALT)  <- gsub("_", ".",
                             colnames(designALT))

fitMALT <- glmFit(GM, designALT) 

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
    levels=designALT)

fitLRT.ALT <- glmLRT(fitMALT, contrast = my.contrasts)

ALL.top.ALT <- topTags(fitLRT.ALT, 1000000)[[1]]

fitLRT.ALT.list <- lapply(colnames(my.contrasts),
                          function (x) glmLRT(fitMALT, contrast = my.contrasts[,x]))
names(fitLRT.ALT.list) <- colnames(my.contrasts)

top.ALT.list <- lapply(fitLRT.ALT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
names(top.ALT.list) <- colnames(my.contrasts)

gene.ALT.list <- lapply(top.ALT.list, function(x) {
                            rownames(x[x$FDR<0.05,])
                        })


