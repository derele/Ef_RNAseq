## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

#setwd("/home/totta/Ef_RNAseq")

if(!exists("raw.counts.4.bg")){
	source("2_edgeR_diff.R")
        }

Ef.RC <- raw.counts.4.bg(Ef.bg)

pdf("figures/rep_pairs_EfG.pdf", width = 70, height=70) # So far not working - because R interrupt?
ggpairs(Ef.RC[[3]])
dev.off()

## use the alternate designE like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
#MOUSE: design <- model.matrix(~0+pData(Ef.bg)$grouped)

## For EIMERIA:
## day 0 rm from pData and designE
designE <- model.matrix(~0+pData(Ef.bg)$grouped)# 0+ in the model formula is an
						# instruction not to include an intercept column and
						# instead to include a column for each group

colnames(designE)  <- gsub("pData\\(Ef.bg\\)\\$grouped", "",
                             colnames(designE))
## just to have the standard point delimitor
colnames(designE)  <- gsub("_", ".",
                             colnames(designE))

keep <- rowSums(Ef.RC[[3]])>100 ## filter lowly expressed
GM.E <- DGEList(Ef.RC[[3]][keep,], group = pData(Ef.bg)$grouped)

GM.E <- calcNormFactors(GM.E)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM.E <- estimateDisp(GM.E, design = designE, robust = TRUE)

pdf("figures/rep_pairs_NORM_EfG.pdf", width = 70, height=70)
ggpairs(cpm(GM.E)) #ggpairs makes matrix of plots
dev.off()

pdf("figures/Ef_genes_mds.pdf")
plotMDS(GM.E, )
dev.off()

pdf("figures/Ef_genes_hclust.pdf")
plot(hclust(dist(t(cpm(GM.E)))))
dev.off()

fitE <- glmFit(GM.E, designE) 

my.contrastsE <- makeContrasts(
    ## just infections in NMRI
#    N3vs0 = NMRI.1stInf.3dpi-NMRI.1stInf.0dpi,
#    N5vs0 = NMRI.1stInf.5dpi-NMRI.1stInf.0dpi,
#    N7vs0 = NMRI.1stInf.7dpi-NMRI.1stInf.0dpi,
    ## just infections in Black and Rag
#    B5vs0 = C57BL6.1stInf.5dpi- C57BL6.0dpi ,
#    R5vs0 = Rag.1stInf.5dpi - Rag.0dpi,
    ## differences in 1st vs 2nd infection
#    N0.1stvsN0.2nd = NMRI.1stInf.0dpi - NMRI.2ndInf.0dpi,
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR2.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,
    ## differences in mouse strains
#        NvsB=NMRI.1stInf.0dpi-C57BL6.0dpi, 
#        NvsR=NMRI.1stInf.0dpi-Rag.0dpi,
#        RvsB=Rag.0dpi-C57BL6.0dpi,
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
    levels=designE)

fitLRT <- glmLRT(fitE, contrast = my.contrastsE)

ALL.top <- topTags(fitLRT, 1000000)

fitLRT.list <- lapply(colnames(my.contrastsE),
                          function (x) glmLRT(fitE, contrast = my.contrastsE[,x]))
names(fitLRT.list) <- colnames(my.contrastsE)

top.list <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
names(top.list) <- colnames(my.contrastsE)

gene.list.E <- lapply(top.list, function(x) {
                            rownames(x[x$FDR<0.01,])
                        })
