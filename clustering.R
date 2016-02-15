## Making heatmaps and other clusters of gene expression data

library(made4)
library(RColorBrewer)

## pheatmap specifying distance and method
pheatmap(cpm(Ef.RC[[3]][unlist(gene.list.E),]),
  color = brewer.pal(n = 8, name = "RdYlBu"),
  show_rownmes = FALSE,
  scale = "row",
  kmeans_k = NA,
 # cluster_rows = T,
 # cluster_cols = T,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  #annotation_names_row = F,
 #cols = brewer.pal(n = 8, name = "PuOr"),
  filename = "figures/pheatmap_Ef_DEG.pdf"
  )
dev.off()

## heatplot is based on heatplot.2 but set to use Pearson correlations as default for distances
pdf("figures/heatplot_Ef.DEG_PuOr.pdf",
  width = 70, #10*300,
  height = 70) #10*300,  
#  res = 300)

heatplot(as.matrix(cpm(Ef.RC[[3]][unlist(gene.list.E), ])),
  scale = "row",
  #margin = c(15, 15), # larger column margin so that labels are not cut off
  labRow = NA,
  cols.default = F,
  #lowcol = "orange", highcol = "purple"
  col = brewer.pal(n = 8, name = "PuOr")
  #annotation_names_row = F,
  )
dev.off()

## compare heatmap.2 which uses Euclidian distances
# Define custom layout for heatmap
mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0)) # creates 3x3 table with location of heatmap elements defined
mylwid = c(1.5,4,0.5)
mylhei = c(1.5,4,1)

pdf("figures/heatmap_Ef.DEG.pdf",
  width = 70,
  height= 70)  
  #res=300)
heatmap(as.matrix(cpm(Ef.RC[[3]][unlist(gene.list.E), ])), 
  scale = "row",
  lmat = mylmat, lwid = mylwid, lhei = mylhei, # to allow key (color legend) to be next to map
  margin = c(15, 15), # larger column margin so that labels are not cut off
  labRow = NULL,
  col = brewer.pal(n = 8, name = "PuOr"),
  keysize = 1.5)
dev.off()


#########################################################
heatmap.2(as.matrix(Ef.RC[[3]][unlist(gene.list.E), ], 
  scale = "row",
  lmat = mylmat, lwid = mylwid, lhei = mylhei, # to allow key (color legend) to be next to map
  #margin = c(5, 1), # larger column margin so that labels are not cut off
  labRow = NA,
  #col = bluered,
  keysize = 1.5))
dev.off()
