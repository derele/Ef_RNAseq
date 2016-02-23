library(RUVSeq)

## We now need a "DGEList" containing edgeR normalized (default TMM)
## data
if(!exists("Mm.RC")){
    source("2_edgeR_diff.R")
}

Mm.differences <- makeGroups(pData(Mm.bg)$group)

## use the expression set and the filtering decided on before
Mm.RUVset <- newSeqExpressionSet(Mm.RC[[3]][keep,])

pdf("figures/ruv_eval/RLE_wo_NORM.pdf")
plotRLE(Mm.RUVset)
dev.off()

pdf("figures/ruv_eval/PCA_wo_NORM.pdf")
plotPCA(Mm.RUVset)
dev.off()

pdf("figures/ruv_eval/mds_wo_NORM.pdf")
plotMDS(counts(Mm.RUVset))
dev.off()

pdf("figures/ruv_eval/hclust_wo_NORM.pdf")
plot(hclust(dist(t(counts(Mm.RUVset)))))
dev.off()

Mm.RUVset.BLnorm <- betweenLaneNormalization(Mm.RUVset, which="upper")

pdf("figures/ruv_eval/RLE_BL_NORM.pdf")
plotRLE(Mm.RUVset.BLnorm)
dev.off()

pdf("figures/ruv_eval/PCA_BL_NORM.pdf")
plotPCA(Mm.RUVset.BLnorm)
dev.off()

pdf("figures/ruv_eval/mds_BL_NORM.pdf")
plotMDS(counts(Mm.RUVset.BLnorm))
dev.off()

pdf("figures/ruv_eval/hclust_BL_NORM.pdf")
plot(hclust(dist(t(normCounts(Mm.RUVset.BLnorm)))))
dev.off()

Mm.RUVset.groups <- RUVs(Mm.RUVset.BLnorm,
                         rownames(Mm.RUVset.BLnorm), k=1, Mm.differences)

pdf("figures/ruv_eval/RLE_group_NORM.pdf") ## not very telling?!
plotRLE(Mm.RUVset.groups)
dev.off()

pdf("figures/ruv_eval/PCA_group_NORM.pdf") ## wow that seems to work a treat
plotPCA(Mm.RUVset.groups)
dev.off()

pdf("figures/ruv_eval/hclust_group_NORM.pdf") ## to much variation removed?
plot(hclust(dist(t(normCounts(Mm.RUVset.groups)))))
dev.off()

pdf("figures/ruv_eval/mds_group_NORM.pdf")
plotMDS(normCounts(Mm.RUVset.groups))
dev.off()

## empirically non-significant genes
Mm.Nsig.genes <- rownames(ALL.top.Mm$table[ALL.top.Mm$table$PValue>0.5, ])

## using empirically ns genes for normalization
Mm.RUVset.empirical <- RUVg(Mm.RUVset.BLnorm, Mm.Nsig.genes, k=1)

pdf("figures/ruv_eval/PCA_empirical_NORM.pdf") ## not as good
plotPCA(Mm.RUVset.empirical)
dev.off()

pdf("figures/ruv_eval/hclust_empirical_NORM.pdf")
plot(hclust(dist(t(normCounts(Mm.RUVset.empirical)))))
dev.off()

pdf("figures/ruv_eval/mds_empirical_NORM.pdf")
plotMDS(normCounts(Mm.RUVset.empirical))
dev.off()

## combination of both
Mm.RUVset.empirical.group <- RUVs(Mm.RUVset.BLnorm, Mm.Nsig.genes, k=1, Mm.differences)

pdf("figures/ruv_eval/PCA_empiricalgroup_NORM.pdf") ## not as good
plotPCA(Mm.RUVset.empirical.group)
dev.off()

pdf("figures/ruv_eval/hclust_empiricalgroup_NORM.pdf")
plot(hclust(dist(t(normCounts(Mm.RUVset.empirical.group)))))
dev.off()

pdf("figures/ruv_eval/mds_empiricalgroup_NORM.pdf")
plotMDS(normCounts(Mm.RUVset.empirical.group))
dev.off()

