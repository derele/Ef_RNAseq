## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

library(limma)
library(scales)
library(stringr)
library(edgeR)

excluded.samples <- c("NMRI_2ndInf_3dpi_rep1",
                      "NMRI_2ndInf_5dpi_rep2",
                      "NMRI_1stInf_0dpi_rep1")

## Mouse and Eimeria specific NORMALIZED counts
Mm.RC <- read.table("output_data/Mm_norm_counts.csv", sep=",")
Ef.RC <- read.table("output_data/Ef_norm_counts.csv", sep=",")

pData <- read.table("output_data/Sample_pData.csv", sep=",")
pData <- pData[!pData$sample%in%excluded.samples, ]

## Mouse and Eimeria specific
Ef.pData <- pData[!pData$dpi%in%"0dpi", ]
Mm.pData <- pData[!pData$challenged%in%c("sporozoites", "oocysts"), ]

## Make sure the pData follows the same order as the sample count
## columns!!!
Mm.RC <- Mm.RC[, order(colnames(Mm.RC))]
Ef.RC <- Ef.RC[, order(colnames(Ef.RC))]

Mm.table <- Mm.pData[order(Mm.pData$sample), ]
Ef.table <- Ef.pData[order(Ef.pData$sample), ]

#################################################################
# MDS CLUSTERING, EIMERIA and MOUSE
######## MOUSE #########################################################
pdf("Supplement/FigureS1_MDS_Mouse.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Mm.RC, labels = Mm.pData$batch, col.axis = "#474747",
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Batch, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$seq.method, col.axis = "#474747",
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Sequencing method, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$dpi,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Dpi, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$challenged, col.axis = "#474747",
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Infection, mouse", line = 0.7)
dev.off()

####### EIMERIA ##########################################################
pdf("Supplement/FigureS1_MDS_Eimeria.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Ef.RC, labels = Ef.pData$batch,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Batch, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$seq.method,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Sequencing method, Ef", line = 0.7)


plotMDS(Ef.RC, labels=Ef.pData$dpi, col.axis = "#474747",
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Dpi, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$challenged,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
title("Infection, Ef", line = 0.7)
dev.off()
