## Making heatmaps and other clusters of gene expression data

library(pheatmap)
library(made4)
library(RColorBrewer)
library(VennDiagram)

## raw count data
Mm.nRC <- read.table("output_data/Mm_norm_counts.csv", sep=",")
Ef.nRC <- read.table("output_data/Ef_norm_counts.csv", sep=",")

## replace zeros and avoid zero variance in cluter scaling
#Mm.nRC[Mm.nRC==0] <- sample(seq(0.01, 0.1, 0.0001), length(Mm.nRC[Mm.nRC==0]), replace =TRUE)
#Ef.nRC[Ef.nRC==0] <- sample(seq(0.01, 0.1, 0.0001), length(Ef.nRC[Ef.nRC==0]), replace =TRUE)

## data on statistical testing for selection
Mm.DE.test <- read.table("output_data/Mm_DEtest.csv", sep=",")
Ef.DE.test <- read.table("output_data/Ef_DEtest.csv", sep=",")

select.diff.data <- function (data, diff.data, columns, FDR.min=0.01){
    DE.test.col <-
        diff.data[diff.data$contrast%in%columns, ]
    DE.diff.genes <-
        unique(DE.test.col[DE.test.col$FDR<FDR.min, "gene"])
    data[DE.diff.genes, ]
}

############# Eimeria LIFE CYCLE SELECTION ##################
cycle.cols <- c("N3vsN5","N3vsN7", "N3vsOoc", "N3vsSpo",
                "N5vsN7", "N5vsOo", "N5vsSp",
                "N7vsOo", "N7vsSp",
                "SpvsOo")
Ef.cycle.diff.data <- select.diff.data(Ef.nRC, Ef.DE.test, cycle.cols, 0.01)

############# Eimeria 1st VERSUS 2nd INFECTION ##################
first.second.cols <- c("B5.1stvsB5.2nd", "N3.1stvsN3.2nd", "N5.1stvsN5.2nd",
                       "N7.1stvsN7.2nd","R5.1stvsR5.2nd")
Ef.first2nd.diff.data <- select.diff.data(Ef.nRC, Ef.DE.test, first.second.cols, 0.5)


############# Eimeria Competent vs. RAG  ##################
comp.rag.cols <- c("R5.1stvsR52ndVSB5.1stvsB52nd", "R5.1stvsR52ndVSN5.1stvsN52nd")
Ef.comp.rag.diff.data <- select.diff.data(Ef.nRC, Ef.DE.test, comp.rag.cols, 0.5)

############# Mouse LIFE CYCLE SELECTION ##################
Mm.cycle.diff.data <- select.diff.data(Mm.nRC, Mm.DE.test, cycle.cols, 0.01)

############# Mouse 1st VERSUS 2nd INFECTION ##################
Mm.first2nd.diff.data <- select.diff.data(Mm.nRC, Mm.DE.test, first.second.cols, 0.01)

############# Mouse Competent vs. RAG  ##################
Mm.comp.rag.diff.data <- select.diff.data(Mm.nRC, Mm.DE.test, comp.rag.cols, 0.2)

get.scaled.and.clustered <- function(data){
    rMeans <- rowMeans(data)
    rSD <- apply(data, 1, sd)
    scaled.data <- scale(t(data), rMeans, rSD)
    hclustered <- hclust(dist(t(scaled.data)))
    return(hclustered)
}

get.cluster.tree.df <- function(data, k){
   clusdat <- get.scaled.and.clustered(data)
   Cluster <- factor(unlist(cutree(clusdat, k = k)))
   as.data.frame(Cluster)
}

Ef.hclustered.cycle <- get.cluster.tree.df(Ef.cycle.diff.data, 7)
Ef.hclustered.first.second <- get.cluster.tree.df(Ef.first2nd.diff.data,  7)

write.table(Ef.hclustered.cycle, "output_data/Ef_hclustered_cycle.csv", sep=",")

Mm.hclustered.cycle <- get.cluster.tree.df(Mm.cycle.diff.data, 7)
Mm.hclustered.first.second <- get.cluster.tree.df(Mm.first2nd.diff.data,  4)
Mm.hclustered.strain <- get.cluster.tree.df(Mm.comp.rag.diff.data, 4)

write.table(Mm.hclustered.cycle, "output_data/Mm_hclustered_cycle.csv", sep=",")

pdf("Supplement/SI_Ef1st2ndHeatmap.pdf",
    height = 8, width = 8, onefile = FALSE) # onefile command to hack away empty page in pdf
pheatmap(Ef.first2nd.diff.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = Ef.hclustered.first.second,
         cutree_rows = 7, 
         show_rownames = F,
         main = expression(paste(italic("E. falciformis"),
             " mRNAs differently abundant between 1st and 2nd infection")))
dev.off()

pdf("figures/Figure3b_EfLifecycleHeatmap.pdf",
    height = 8, width = 8, onefile = FALSE) # onefile to rm empty page in pdf
pheatmap(Ef.cycle.diff.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = Ef.hclustered.cycle,
         cutree_rows = 7, 
         show_rownames = F,
         main = expression(paste(italic("E. falciformis"),
                                 " mRNAs differently abundant between lifecycle stages"))
         )
dev.off()

##############################################
## Mouse plotting
pdf("Supplement/SI_MmStrainHeatmap.pdf",
    height = 8, width = 8, onefile = FALSE) # onefile to rm empty page in pdf
pheatmap(Mm.comp.rag.diff.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = Mm.hclustered.strain,
         ## annotation_names_row = F,
         cutree_rows = 4, 
         show_rownames = F,
         main = expression(paste(italic("M. musculus"),
             " mRNAs differently abundant in different mouse strains")))
dev.off()

pdf("Supplement/SI_Mm1st2ndHeatmap.pdf",
    height = 8, width = 8, onefile = FALSE) # onefile to rm empty page in pdf
pheatmap(Mm.first2nd.diff.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = Mm.hclustered.first.second,
         cutree_rows = 4, 
         show_rownames = F,
         main = expression(paste(italic("M. musculus"),
             " mRNAs differently abundant in 1st and2nd infection")))
dev.off()


pdf("figures/Figure2b_MmLifecycleHeatmap.pdf",
    height = 8, width = 8, onefile = FALSE) # onefile to rm empty page in pdf
pheatmap(Mm.cycle.diff.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         annotation_row = Mm.hclustered.cycle,
         cutree_rows = 7, 
         show_rownames = F,
         main = expression(paste(italic("M. musculus"),
             " mRNAs differently abundant at different dpi ")))
dev.off()

###########################
## Check for specific genes withing clusters
efab.clus1 <- subset(Ef.hclustered.cycle, Cluster=='1')
efab.clus2 <- subset(Ef.hclustered.cycle, Cluster=='2')
efab.clus3 <- subset(Ef.hclustered.cycle, Cluster=='3')
efab.clus4 <- subset(Ef.hclustered.cycle, Cluster=='4')
efab.clus5 <- subset(Ef.hclustered.cycle, Cluster=='5')
efab.clus6 <- subset(Ef.hclustered.cycle, Cluster=='6')
efab.clus7 <- subset(Ef.hclustered.cycle, Cluster=='7')

write.table(efab.clus1, "Supplement/SI_genes_efCluster1.csv", sep=",")
write.table(efab.clus2, "Supplement/SI_genes_efCluster2.csv", sep=",")
write.table(efab.clus3, "Supplement/SI_genes_efCluster3.csv", sep=",")
write.table(efab.clus4, "Supplement/SI_genes_efCluster4.csv", sep=",")
write.table(efab.clus5, "Supplement/SI_genes_efCluster5.csv", sep=",")
write.table(efab.clus6, "Supplement/SI_genes_efCluster6.csv", sep=",")
write.table(efab.clus7, "Supplement/SI_genes_efCluster7.csv", sep=",")

####################################
## check for surface antigen genes in Eimeria clusters (SAG-genes)
####################################

## gene IDs
sag.genes <- c("EfaB_PLUS_193.g6",
               "EfaB_PLUS_49504.g2703",
               "EfaB_PLUS_25458.g2090",
               "EfaB_PLUS_17089.g1451",
               "EfaB_PLUS_3222.g332",
               "EfaB_MINUS_7048.g636",
               "EfaB_MINUS_56725.g2963",
               "EfaB_PLUS_7048.g677",
               "EfaB_PLUS_2074.g242",
               "EfaB_MINUS_22450.g1915",
               "EfaB_MINUS_7048.g637",
               "EfaB_MINUS_6743.g610")

sag.check <- subset(efab.clus1, row.names(efab.clus1) %in% sag.genes)
sag.check <- rbind(sag.check, subset(efab.clus2, row.names(efab.clus2) %in% sag.genes))
sag.check <- rbind(sag.check, subset(efab.clus3, row.names(efab.clus3) %in% sag.genes))
sag.check <- rbind(sag.check, subset(efab.clus4, row.names(efab.clus4) %in% sag.genes))
sag.check <- rbind(sag.check, subset(efab.clus5, row.names(efab.clus5) %in% sag.genes))
sag.check <- rbind(sag.check, subset(efab.clus6, row.names(efab.clus6) %in% sag.genes))
sag.check <- rbind(sag.check, subset(efab.clus7, row.names(efab.clus7) %in% sag.genes))

################### ROPs #####################

rop.genes <- c("EfaB_MINUS_17096.g1521",	
               "EfaB_MINUS_32658.g2475",	
               "EfaB_MINUS_42996.g2710",
               "EfaB_MINUS_720.g57",
               "EfaB_PLUS_15899.g1411",		
               "EfaB_PLUS_24117.g1969",
               "EfaB_PLUS_33184.g2393",
               "EfaB_PLUS_47595.g2679",
               "EfaB_PLUS_7742.g778",
               "EfaB_PLUS_8664.g829")

rop.check <- subset(efab.clus1, row.names(efab.clus1) %in% rop.genes)
rop.check <- rbind(rop.check, subset(efab.clus2, row.names(efab.clus2) %in% rop.genes))
rop.check <- rbind(rop.check, subset(efab.clus3, row.names(efab.clus3) %in% rop.genes))
rop.check <- rbind(rop.check, subset(efab.clus4, row.names(efab.clus4) %in% rop.genes))
rop.check <- rbind(rop.check, subset(efab.clus5, row.names(efab.clus5) %in% rop.genes))
rop.check <- rbind(rop.check, subset(efab.clus6, row.names(efab.clus6) %in% rop.genes))
rop.check <- rbind(rop.check, subset(efab.clus7, row.names(efab.clus7) %in% rop.genes))


