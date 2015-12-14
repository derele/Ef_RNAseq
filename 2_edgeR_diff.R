#### DIFFERENTIAL GENE EXPRESSION ANALYSIS
#### OF RNA-seq DATA
## LIBRARIES
library(edgeR)
library(RSvgDevice)

raw.counts.4.bg <- function(bg){
exon.raw.count <- eexpr(bg, "ucount")
## the linkage data
e2t <- bg@indexes$e2t
t2g <- bg@indexes$t2g 
e2t2g <- merge(e2t, t2g)
all.count <- merge(e2t2g, exon.raw.count,
                   by.x = "e_id", by.y = 0)
## sum up for transcripts
transcript.raw.count <-
    do.call("rbind",
            by(all.count, all.count$t_id,
               function(x) colSums(x[, 4:ncol(x)])))
## sum up for genes
gene.raw.count <-
    do.call("rbind",
            by(all.count, all.count$g_id,
               function(x) colSums(x[, 4:ncol(x)])))
colnames(transcript.raw.count) <-
    as.character(pData(bg)$samples)
colnames(gene.raw.count) <-
    as.character(pData(bg)$samples)
return(list(exon.raw.count,
            transcript.raw.count,
            gene.raw.count))
}

mouse.RC <- raw.counts.4.bg(mouse.bg)

design <- model.matrix(~pData(mouse.bg)$infection +
                           pData(mouse.bg)$mouse.strain +
                               pData(mouse.bg)$timepoint)

colnames(design)  <- gsub("pData\\(mouse.bg\\)\\$", "",
                          colnames(design))

GM <- DGEList(mouse.RC[[3]], group = pData(mouse.bg)$grouped)
GM <- calcNormFactors(GM)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM <- estimateDisp(GM, design = design, robust = TRUE)

## devSVG("figures/Mouse_genes_mds.svg")
plotMDS(GM)
## dev.off()

## devSVG("figures/Mouse_genes_hclust.svg")
## plot(hclust(dist(t(cpm(GM)))))
## dev.off()

fitM <- glmFit(GM, design) 
fitLRT.list <- lapply(1:6, function (i) glmLRT(fitM, coef=i))

list.of.top.genesM <- lapply(fitLRT.list, topTags, 1000000)

names(list.of.top.genesM) <- unlist(lapply(list.of.top.genesM, 
                                    function(x) x@.Data[[3]]))

list.of.diff.genes <- lapply(list.of.top.genesM, function(x) {
                                 d  <- x[[1]]
                                 rownames(d[d$FDR<0.05,])
                             })


length(list.of.diff.genes[[2]])
length(chal_genes)
length(chal_genesI)

length(intersect(list.of.diff.genes[[2]], chal_genes))
length(intersect(list.of.diff.genes[[2]], chal_genesI))






