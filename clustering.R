## Making heatmaps and other clusters of gene expression data

library(made4)
library(RColorBrewer)


########## USE THIS #######################
########## PHEATMAP() specifying distance and method ############
pheatmap(cpm(Ef.RC[[3]][unlist(gene.list.E),]),
  color = brewer.pal(n = 8, name = "RdBu"),
  scale = "row",
 # kmeans_k = NA,
 # cluster_rows = T,
 # cluster_cols = T,
 # clustering_distance_rows = "euclidean", --> has no influence on clustering
 # clustering_method = "complete",         --> has no influence on clustering
 #annotation_names_row = F,
  annotation_row = NA,
  show_rownames = F,
  main = "E. falciformis, (complete, euclidean)",
 #cols = brewer.pal(n = 8, name = "PuOr"),
  filename = "figures/pheatmap_Ef_DEGy.pdf"
  )
dev.off()

########### MAKE SUBSET OF DEG dataset for testing heatmap functions #############
#set.seed(111) #to get same output from sample
#efsample <- sapply(c(10:15, 20:25), function (e) rnbinom(20, size=e, prob = 0.05))
#rownames(efsample) <- 1:nrow(efsample)
#colnames(efsample) <- 1:ncol(efsample)

#pheatmap(efsample)
#heatmap(efsample)
### CONCLUSION:
# heatmap calls as.dendrogram() and reorderfun() (internal function?) before plotting, whereas pheatmap plots the
# hclust object without it. 
# as.dendrogram() does some kind of sorting (or is it reorderfun() only?) - find out what!

### OBJECT TO USE FOR HEATMAPS AND PHEATMAPS, but doesn't work (yet?)
#Efdist <- dist(efsample) #[unlist(gene.list.E),]))
#Efhclust <- hclust(Efdist) # because heatmap functions require numberic matrix
#Efdend <- as.dendrogram(Efhclust)
#plot(Efdend, horiz=TRUE)

#pheatmap(cpm(Ef.RC[[3]][unlist(gene.list.E),]), clustering_distance_rows = Efdist, clustering_distance_cols = t(Efdist))

################## HEATMAP() as comparison ##################
#pdf("figures/heatmap_Ef_DEGa.pdf",
#  width = 10,
#  height= 10)  
## heatmap() uses hclust (default 'complete') and dist (default 'euclidean')
#heatmap(cpm(Ef.RC[[3]][unlist(gene.list.E), ]), 
#  hclustfun = Efdist
#  Rowv = NA, Colv = NA,
#  scale = "row",
#  margin = c(15, 5),
#  #cexRow = 10,
#  labRow = NA, # NA to hide rownames
#  main = "heatmap(), Ef, 'complete', 'euclidean'",
#  col = brewer.pal(n = 8, name = "RdBu"))
#dev.off()

# Define custom layout for heatmap - can probably be removed
#mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0)) # creates 3x3 table with location of heatmap elements defined
#mylwid = c(1.5,4,0.5)
#mylhei = c(1.5,4,1)
# lmat = mylmat, lwid = mylwid, lhei = mylhei, # to allow key (color legend) to be next to map

