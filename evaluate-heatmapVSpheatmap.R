## testing if heatmap() and pheatmap() really cluster the same way
## initially looked funky, but now I think there is no difference
## some forum comments also indicated that the callback-object below
## is needed for pheatmap to reproduce heatmap,
## but in the end I could not see this difference either, i.e. they
## cluster identically.

library(made4)
library(RColorBrewer)

## dummy neg binom data
efsample <- sapply(c(10:15, 20:25), function (e) rnbinom(20, size=e, prob = 0.05))
rownames(efsample) <- 1:nrow(efsample)
colnames(efsample) <- 1:ncol(efsample)

##### function from pheatmap manual to set callback to reproduce heatmap() row sorting
callback = function(hc, mat){
sv = svd(t(mat))$v[,1]
dend = reorder(as.dendrogram(hc), wts = sv)
as.hclust(dend)
}
## pheatmap test --> run with and without cluster_callback spec
pheatmap(efsample, 
        #cluster_callback = callback,
        color = brewer.pal(n = 8, name = "RdBu"),
        main = "pheatmap test, default settings",
        filename = "figures/pheat-test-default.pdf"
        )
dev.off()
#### heatmap test
pdf("figures/heat-test.pdf", 
    width=10, height=10)
heatmap(efsample,
        col = brewer.pal(n = 8, name = "RdBu"),
        main = "heatmap test, default"
        )
dev.off()

### CONCLUSION:
# heatmap calls as.dendrogram() and reorderfun() (internal function?) before plotting, whereas pheatmap plots the
# hclust object without it. 
# There is however no meaningful (interpretation-wise) difference between the two.


