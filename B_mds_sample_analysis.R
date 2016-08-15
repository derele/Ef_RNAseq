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
pData <- pData[!pData$sample%in%excluded.samples,]

## Mouse and Eimeria specific
Ef.pData <- pData[!pData$dpi%in%0, ]

## Make sure the pData follows the same order as the sample count
## columns!!!
Mm.RC <- Mm.RC[, order(colnames(Mm.RC))]
Ef.RC <- Ef.RC[, order(colnames(Ef.RC))]

Mm.pData <- Mm.pData[order(Mm.pData$sample), ]
Ef.pData <- Ef.pData[order(Ef.pData$sample), ]

###################################################################
## COLORS FOR pDATA() GROUPS 
###################################################################

batch.colors.ef = c("#1b9e77", # 3
		 "#1b9e77", # 3
		 "#1b9e77", # 3
                  "#d95f02", # 1
                  "#003366", # 100
                  "#9150f9", # 2
                  "#9150f9", # 2
                  "#474747", # 0
		  "#d95f02", # 1
		  "#d95f02", # 1
                  "#1b9e77", # 3
                  "#9150f9", # 2
                  "#1b9e77",# 3
		  "#1b9e77", # 3
                  "#d95f02", # 1
                  "#474747", # 0
		  "#d95f02", # 1
                  "#474747", # 0
                  "#1b9e77", "#1b9e77", "#1b9e77") # 3

seq.colors.ef = c("#474747", "#474747", "#474747",
               "#1b9e77", "#1b9e77", "#1b9e77",
               "#1b9e77","#1b9e77", "#1b9e77", "#1b9e77",
               "#474747",
               "#1b9e77",
               "#474747","#474747",
               "#1b9e77", "#1b9e77", "#1b9e77", "#1b9e77",
               "#474747","#474747","#474747")

batch.colors.mm = c("#1b9e77", # 3
		 "#1b9e77", # 3
		 "#1b9e77", # 3
		 "#1b9e77", # 3
		 "#1b9e77", # 3
                  "#d95f02", # 1
                  "#003366", # 100
                  "#9150f9", # 2
                  "#9150f9", # 2
                  "#474747", # 0
                  "#d95f02", # 1
                  "#d95f02", # 1
                  "#9150f9", # 2
		  "#1b9e77", # 3
		  "#1b9e77", # 3
                  "#9150f9", # 2
		  "#1b9e77", # 3
		  "#1b9e77", # 3
                  "#1b9e77", "#1b9e77", "#1b9e77","#1b9e77", "#1b9e77") # 3

seq.colors.mm = c("#474747", "#474747", "#474747", "#474747", "#474747",
               "#1b9e77", "#1b9e77", "#1b9e77","#1b9e77",
               "#1b9e77","#1b9e77", "#1b9e77", "#1b9e77",
               "#474747","#474747",
               "#1b9e77",
               "#474747","#474747","#474747", "#474747",
               "#474747","#474747","#474747")

#################################################################
# MDS CLUSTERING, EIMERIA and MOUSE
#################################################################
pdf("Supplement/MDS_Mouse.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Mm.RC, labels = Mm.pData$batch, col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = batch.colors.mm,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$seq.method, col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = seq.colors.mm,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequencing method, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$dpi,
        ## col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        col = batch.colors.mm,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Dpi, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$challenged, col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = seq.colors.mm,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Infection, mouse", line = 0.7)
dev.off()



#################################################################
pdf("Supplement/MDS_Eimeria.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Ef.RC, labels=Ef.pData$batch,
        ## col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = batch.colors.ef,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$seq.method,
        ## col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = seq.colors.ef,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequencing method, Ef", line = 0.7)


plotMDS(Ef.RC, labels=Ef.pData$dpi, col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = batch.colors.ef,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Dpi, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$challenged,
        ## col.axis = "#474747",
        ## col.lab = "#474747",
        ## col.main = "#474747",
        ## col.sub = "#474747",
        ## col = seq.colors.ef,
        xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Infection, Ef", line = 0.7)


dev.off()
