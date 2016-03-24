## Making heatmaps and other clusters of gene expression data

library(pheatmap)
library(made4)
library(RColorBrewer)


########## USE THIS #######################
########## PHEATMAP() specifying distance and method ############

union.of.Ef.cycle.diff.top.100 <-
    unique(unlist(lapply(Ef.1st.pass.model[[3]][1:10], head, n=500)))

Ef.cycle.diff.top.100.data <-
    cpm(Ef.1st.pass.model[[4]])[union.of.Ef.cycle.diff.top.100,]

get.scaled.and.clustered <- function(data){
    rMeans <- rowMeans(data)
    rSD <- apply(data, 1, sd)
    scaled.data <- scale(t(data), rMeans, rSD)
    hclustered <- hclust(dist(t(scaled.data)))
    return(hclustered)
}

Ef.hclustered <- get.scaled.and.clustered(Ef.cycle.diff.top.100.data)

hcluster <- list()
hcluster[["Ef"]] <- cutree(Ef.hclustered, k= 5) ## h=7.2)

Ef.hclustered.df <- as.data.frame(hcluster[["Ef"]])
names(Ef.hclustered.df) <- "Cluster"
Ef.hclustered.df$Cluster <- as.factor(Ef.hclustered.df$Cluster)

pdf("/home/ele/Dropbox/Totta_tmp_Efposter/Ef_most_sig_lifecycle_heatmap.pdf",
    height = 8, width = 8)
pheatmap(Ef.cycle.diff.top.100.data,
         color = brewer.pal(n = 11, name = "BrBG"), 
         scale = "row",
         cluster_rows = T, ## hc.high,
         cluster_cols = T,
         annotation_row = Ef.hclustered.df,
         ## annotation_names_row = F,
         cutree_rows = 5, 
  show_rownames = F,
         main = expression(paste(italic("E. falciformis"),
             " genes differentially expressed between lifecycle stages")))
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

pdf("/home/ele/Dropbox/Totta_tmp_Efposter/Mm_most_sig_lifecycle_heatmap.pdf",
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
             " genes differentially expressed at different dpi ")))
dev.off()

