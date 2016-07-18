library(ggplot2)
library(reshape)
library(plyr)
library(phangorn)
#library(multicore)

if(!exists("omcl")){
  source("3_annotations.R")
}

if(!exists("Ef.1st.pass.model")){
  source("2_edgeR_diff.R")
}


## new 2016: Get for every cluster the E. falciformis and E. tenella gene
Efa.Ete <- lapply(omcl, function (x) {
                      x[grep("Efa\\||Ete\\|", x)]})

## has to be a cluster still 
Efa.Ete <- Efa.Ete[unlist(lapply(Efa.Ete, length))>1]

## and contain genes from both
Efa.Ete <- Efa.Ete[unlist(lapply(Efa.Ete, function (x){ any(grepl("Efa\\|", x)) &
                                                           any(grepl("Efa\\|", x)) }))]

## split them appart again
Efa.list <-lapply(Efa.Ete, function (x) x[grep("Efa\\|", x)])
Ete.list <- lapply(Efa.Ete, function (x) x[grep("Ete\\|", x)])

## Do all combinations of each Efa gene with with each Ete gene in a cluster
cobn.list <- lapply(1:length(Efa.list), function(i){
                        exp <- rep(Efa.list[[i]], times=length(Ete.list[[i]]))
                        cbind(exp, Ete.list[[i]])
                    })

all.orthologs <- as.data.frame(do.call("rbind", cobn.list))
names(all.orthologs) <- c("Efa", "Eten")
all.orthologs$Eten <- gsub("Ete\\|", "", all.orthologs$Eten)
all.orthologs$Efa <- gsub("Efa\\|", "", all.orthologs$Efa)


## read the Tenella data
TenGam <- read.csv("data/12864_2015_1298_MOESM3_ESM.csv", skip=1)

names(TenGam)[5:8] <- paste(names(TenGam)[5:8], "GamVsMer", sep=".")
names(TenGam) <- gsub("\\.1", "\\.GamVsSpor", names(TenGam))

## make thins numeric
TenGam[, c("gametocyte", "merozoite", "sporozoite",
  "fold.change.GamVsSpor", "log2.fold.change.GamVsSpor",
           "fold.change.GamVsMer", "log2.fold.change.GamVsMer")] <-
  apply(TenGam[, c("gametocyte", "merozoite", "sporozoite",
             "fold.change.GamVsSpor", "log2.fold.change.GamVsSpor",
                   "fold.change.GamVsMer", "log2.fold.change.GamVsMer")], 2,
        function (x) as.numeric(as.character(x)))

## Work with normalized counts
RC.ortho <- merge(all.orthologs,
                  TenGam[, c("gene.ID", "gametocyte", "merozoite", "sporozoite")],
                  by.x="Eten", by.y="gene.ID")

RC.ortho <- merge(RC.ortho, cpm(Ef.1st.pass.model[[4]]),
                  by.x="Efa", by.y=0)
RC.ortho[,3:ncol(RC.orhto)] <- apply(RC.ortho[,3:ncol(RC.orhto)], 2,
                                     function (x) as.numeric(as.character(x)))

CS <- colSums(RC.ortho[,3:ncol(RC.ortho)])
scaling.f <- max(CS)/CS

RC.ortho[, 3:ncol(RC.ortho)] <- t(t(RC.ortho[, 3:ncol(RC.ortho)])*scaling.f)

pdf("figuresANDmanuscript/RC_correlation_Efal_Eten.pdf")
pheatmap(cor(RC.ortho[,3:ncol(RC.ortho)], method="spearman")[-(1:3),1:3],
         display_numbers=TRUE)
dev.off()


###
mean.columns <- function(x){
  reps <- as.factor(gsub("_rep\\d", "", names(x)[6:ncol(x)]))
  dat <- t(x[,6:ncol(x)])
  bar <- do.call(rbind, by(dat, reps, colMeans))
  t(bar)
}
  
RC.ortho.mean <- as.data.frame(cbind(RC.ortho[, 1:5], mean.columns(RC.ortho)))

pdf("figuresANDmanuscript/RC_mean_correlation_Efal_Eten.pdf")
pheatmap(cor(RC.ortho.mean[,3:ncol(RC.ortho.mean)], method="spearman")[-(1:3),1:3],
         display_numbers=TRUE)
dev.off()

### Work with Fold Changes... ## What to do wiht Max and Min?
FC.ortho <- merge(all.orthologs,
                  TenGam[, grepl("gene.ID|log2.fold.change", names(TenGam))],
                  by.x="Eten", by.y="gene.ID")
FC.ortho <- FC.ortho[!rowSums(is.na(FC.ortho))>0, ]

all.Ef.FC.list <- lapply(Ef.1st.pass.model[[2]], function(x) {
  as.data.frame(cbind(gene=rownames (x), logFC=x[,1]))
})

all.FC.df <- Reduce(function (x,y)merge(x,y, by="gene"),all.Ef.FC.list)
names(all.FC.df) <- c("gene.ID", paste("logFC", names(Ef.1st.pass.model[[2]]), sep="."))

FC.ortho <- merge(FC.ortho,
                  all.FC.df,
                  by.x="Efa", by.y="gene.ID")

FC.ortho[,3:ncol(FC.ortho)] <- apply(FC.ortho[,3:ncol(FC.ortho)],
                                     2, function(x) as.numeric(as.character(x)))


pdf("figuresANDmanuscript/FC_correlation_Efal_Eten.pdf")
pheatmap(cor(FC.ortho[, 3:ncol(FC.ortho)], method="spearman")[-c(1:2),1:2],
         display_numbers=TRUE)
dev.off()



pdf("figuresANDmanuscript/TenellaGamVsMer_vs_N3vsN7.pdf")
ggplot(FC.ortho, aes(log2.fold.change.GamVsMer, logFC.N3vsN7*-1)) +
  geom_point(alpha=0.8) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
  stat_smooth() +
  theme_bw()
dev.off()


pdf("figuresANDmanuscript/TenellaGamVsMer_vs_N5vsN7.pdf")
ggplot(FC.ortho, aes(log2.fold.change.GamVsMer, logFC.N5vsN7*-1)) +
  geom_point(alpha=0.8) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
  stat_smooth() +
  theme_bw()
dev.off()


pdf("figuresANDmanuscript/TenellaGamVsSpor_vs_N7vsSpor.pdf")
ggplot(FC.ortho, aes(log2.fold.change.GamVsSpor, logFC.N7vsSp)) +
  geom_point(alpha=0.8) +
  stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
  stat_smooth() +
  theme_bw()
dev.off()

