## Making heatmaps and other clusters of gene expression data

library(pheatmap)
library(made4)
library(RColorBrewer)

if(!exists("RNAseq.Array.logFC")){
    source("4_compareArray.R")
}
########## USE THIS #######################
########## PHEATMAP() specifying distance and method ############
## 1:10 refers to the Ef.contrast object and the first ten 
## comparisons in it.
############# LIFE CYCLE SELECTION ##################
union.of.Ef.cycle.diff.top.100 <-
    unique(unlist(lapply(Ef.1st.pass.model[[3]][1:10], head, n=500)))

Ef.cycle.diff.top.100.data <-
    cpm(Ef.1st.pass.model[[4]])[union.of.Ef.cycle.diff.top.100,]

############# 1st VERSUS 2nd INFECTION ##################
union.of.Ef.1st2nd.top <-
    unique(unlist(lapply(Ef.1st.pass.model[[3]][11:15], head, n=500)))

Ef.1st2nd.top.data <-
    cpm(Ef.1st.pass.model[[4]])[union.of.1st2nd.top,]

######## EIMERIA DEVELOPMENT IN DIFFERENT MOUSE STRAINS (day 5 only) ##################
union.of.Ef.strain.depend.top <-
    unique(unlist(lapply(Ef.1st.pass.model[[3]][16:18], head, n=500)))

Ef.strain.depend.top.data <-
    cpm(Ef.1st.pass.model[[4]])[union.of.Ef.strain.depend.top,]

get.scaled.and.clustered <- function(data){
    rMeans <- rowMeans(data)
    rSD <- apply(data, 1, sd)
    scaled.data <- scale(t(data), rMeans, rSD)
    hclustered <- hclust(dist(t(scaled.data)))
    return(hclustered)
}

Ef.hclustered <- get.scaled.and.clustered(Ef.cycle.diff.top.100.data)
Ef.hclustered.first.second <- get.scaled.and.clustered(Ef.1st2nd.top.data) 
Ef.hclustered.strain <- get.scaled.and.clustered(Ef.strain.depend.top.data )


hcluster <- list()
hcluster[["Ef"]] <- cutree(Ef.hclustered, k= 7) ## TK changed k from 5 ## h=7.2)

Ef.hclustered.df <- as.data.frame(hcluster[["Ef"]])
names(Ef.hclustered.df) <- "Cluster"
Ef.hclustered.df$Cluster <- as.factor(Ef.hclustered.df$Cluster)

pdf("figures/Ef_most_sig_lifecycle_heatmap.pdf",
    height = 8, width = 8)
pheatmap(Ef.cycle.diff.top.100.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T, ## hc.high,
         cluster_cols = T,
         annotation_row = Ef.hclustered.df,
         ## annotation_names_row = F,
         cutree_rows = 7, # 5 in EH script 
         show_rownames = F,
         main = expression(paste(italic("E. falciformis"),
             " mRNAs differently abundant between lifecycle stages")))
dev.off()


## mouse 
union.of.Mm.cycle.diff.top.100 <-
    unique(unlist(lapply(Mm.1st.pass.model[[3]][c(1:3, 12:14)], head, n=500)))

Mm.cycle.diff.top.100.data <-
    cpm(Mm.1st.pass.model[[4]])[union.of.Mm.cycle.diff.top.100,]

Mm.hclustered <- get.scaled.and.clustered(Mm.cycle.diff.top.100.data)

hcluster[["Mm"]] <- cutree(Mm.hclustered, k= 4) ## h=7.2)

Mm.hclustered.df <- as.data.frame(hcluster[["Mm"]])

names(Mm.hclustered.df) <- "Cluster"
Mm.hclustered.df$Cluster <- as.factor(Mm.hclustered.df$Cluster)

pdf("figures/Mm_most_sig_lifecycle_heatmap.pdf",
    height = 8, width = 8)
pheatmap(Mm.cycle.diff.top.100.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T, ## hc.high,
         cluster_cols = T,
         annotation_row = Mm.hclustered.df,
         ## annotation_names_row = F,
         cutree_rows = 4, 
         show_rownames = F,
         main = expression(paste(italic("M. musculus"),
             " mRNAs differently abundant at different dpi ")))
dev.off()

