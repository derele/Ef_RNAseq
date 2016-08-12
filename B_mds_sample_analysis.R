## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

library(limma)
library(scales)
library(stringr)

library(edgeR)

excluded.samples <- c("NMRI_2ndInf_3dpi_rep1",
                      "NMRI_2ndInf_5dpi_rep2",
                      "NMRI_1stInf_0dpi_rep1")

All.RC <- as.matrix(read.table("output_data/RC_All_genes.csv", sep=","))

All.RC <- All.RC[, !colnames(All.RC)%in%excluded.samples]

## Mouse and Eimeria specific
Mm.RC <- All.RC[grepl('^ENSMUS.*', rownames(All.RC)),]
Ef.RC <- All.RC[grepl('^EfaB.*', rownames(All.RC)),]

pData <- read.table("output_data/Sample_pData.csv", sep=",")
pData <- pData[!pData$samples%in%excluded.samples, ]

## Mouse and Eimeria specific
Ef.pData <- pData[!pData$timepoint%in%0, ]
Mm.pData <- pData[!is.na(pData$timepoint), ]

Ef.pData$grouped <- as.factor(as.character(Ef.pData$grouped))
Mm.pData$grouped <- as.factor(as.character(Mm.pData$grouped))

## Mouse and Eimeri spcific samples
Mm.RC <- Mm.RC[, colnames(Mm.RC)%in%Mm.pData$samples]
Ef.RC <- Ef.RC[, colnames(Ef.RC)%in%Ef.pData$samples] 

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
pdf("Supplement/MDS_technical.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Ef.RC, labels=Ef.pData$batch, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = batch.colors.ef, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$seq.method,
        col.axis = "#474747", col.lab = "#474747",
        col.main = "#474747", col.sub = "#474747",
        col = seq.colors.ef, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequencing method, Ef", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$batch, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = batch.colors.mm, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$seq.method, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = seq.colors.mm, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequencing method, mouse", line = 0.7)

dev.off()



#################################################################
pdf("figures/figure3_MDS_biology.pdf")
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))
plotMDS(Ef.RC, labels=Ef.pData$timepoint, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = batch.colors.ef, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Dpi, Ef", line = 0.7)

plotMDS(Ef.RC, labels = Ef.pData$infection,
        col.axis = "#474747", col.lab = "#474747",
        col.main = "#474747", col.sub = "#474747",
        col = seq.colors.ef, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Infection, Ef", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$timepoint, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = batch.colors.mm, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Dpi, mouse", line = 0.7)

plotMDS(Mm.RC, labels = Mm.pData$infection, col.axis = "#474747",
        col.lab = "#474747", col.main = "#474747",
        col.sub = "#474747",
        col = seq.colors.mm, xlab = "Fold change, dimension 1",
        ylab = "Fold change, dimension 2", mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Infection, mouse", line = 0.7)

dev.off()
