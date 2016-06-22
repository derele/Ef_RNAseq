## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

#setwd("/home/totta/Ef_RNAseq")

if(!exists("raw.counts.4.bg")){
	source("2_edgeR_diff.R")
        }

library(limma)
library(scales)
library(stringr)

#Ef.RC <- raw.counts.4.bg(Ef.bg)

######################################

######################################
pdf("figures/rep_pairs_EfG.pdf", width = 70, height=70) # So far not working - because R interrupt?
ggpairs(Ef.RC[[3]])
dev.off()


##################### For EIMERIA: #################################

ggplot()+
    geom_point(data = tableEf, aes(x = Samples,
                   y = as.numeric(colSums(Ef.RC[[3]])))) +
	scale_y_log10(labels = comma) +
            ylab("Number of reads") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))



pdf("figures/rep_pairs_NORM_EfG.pdf", width = 70, height=70)
ggpairs(cpm(GM.E)) #ggpairs makes matrix of plots
dev.off()

###################################################################
## COLORS FOR pDATA() GROUPS 
###################################################################

### Show all the RColorBrewer colour schemes available
display.brewer.all()
####

#day.colors = c("#1b9e77", "#1b9e77", "#1b9e77", 
#               "#d95f02", "#d95f02", 
#               "#1b9e77", "#1b9e77", 
#               "#474747", "#474747",
#               "#d95f02", "#d95f02", 
#               "#1b9e77", "#1b9e77",
#               "#474747", "#474747",
#               "#1b9e77","#1b9e77", "#1b9e77")

#strain.colors = c("#1b9e77", "#1b9e77", "#1b9e77",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#474747", "#474747",
#                "#d95f02", "#d95f02", "#d95f02") 
## Batches 0, 1, 2 are non-hiseq and 3 are HU_SS ones. 100 are not clear.
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
png("figures/EfMm_4-mds.png")
#plotMDS(GM.E, )
par(mfrow = c(2,2), mai = c(0.6, 0.5, 0.4, 0.4))

plotMDS(Ef.RC[[3]], labels = pData(Ef.bg)$batch,
	      col.axis = "#474747", col.lab = "#474747", col.main = "#474747", col.sub = "#474747",
              col = batch.colors.ef,
              xlab = "Fold change, dimension 1", ylab = "Fold change, dimension 2",
              mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, Ef", line = 0.7)

plotMDS(Ef.RC[[3]], labels = pData(Ef.bg)$seq.method,          
              col.axis = "#474747", col.lab = "#474747", col.main = "#474747", col.sub = "#474747",
              col = seq.colors.ef,
              xlab = "Fold change, dimension 1", ylab = "Fold change, dimension 2",
              mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequnecing method, Ef", line = 0.7)


plotMDS(Mm.RC[[3]], labels = pData(Mm.bg)$batch,
              col.axis = "#474747", col.lab = "#474747", col.main = "#474747", col.sub = "#474747",
              col = batch.colors.mm,
              xlab = "Fold change, dimension 1", ylab = "Fold change, dimension 2",
              mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1)
title("Batch, mouse", line = 0.7)

plotMDS(Mm.RC[[3]], labels = pData(Mm.bg)$seq.method,          
              col.axis = "#474747", col.lab = "#474747", col.main = "#474747", col.sub = "#474747",
              col = seq.colors.mm,
              xlab = "Fold change, dimension 1", ylab = "Fold change, dimension 2",
              mgp = c(2,1,0))
#mai = c(1, 0.1, 0.1, 0.1))
title("Sequnecing method, mouse", line = 0.7)

dev.off()



