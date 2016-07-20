## Description of script

library(ggplot2)
=======
library(pheatmap)
library(reshape)
library(ggplot2)

if(!exists("omcl")){
  source("3_annotations.R")
}

if(!exists("Ef.1st.pass.model")){
  source("2_edgeR_diff.R")
}


## new 2016: Get for every cluster the E. falciformis and E. tenella gene
get.all.pairedE.orthos <- function (){
    Efa.Ete <- lapply(omcl, function (x) {
        x[grep("Efa\\||Ete\\|", x)]
    })
    ##  to be a cluster still 
    Efa.Ete <- Efa.Ete[unlist(lapply(Efa.Ete, length))>1]
    ## and contain genes from both
    Efa.Ete <- Efa.Ete[unlist(lapply(Efa.Ete, function (x){ any(grepl("Efa\\|", x)) &
                                                                any(grepl("Ete\\|", x)) }))]
    ## split them appart again
    Efa.list <-lapply(Efa.Ete, function (x) x[grep("Efa\\|", x)])
    Ete.list <- lapply(Efa.Ete, function (x) x[grep("Ete\\|", x)])
    ## Do all combinations of each Efa gene with with each Ete gene in a cluster
    cobn.list <- lapply(1:length(Efa.list), function(i){
        exp <- rep(Efa.list[[i]], times=length(Ete.list[[i]]))
        cbind(exp, Ete.list[[i]])
    })
    all.orthologs <- as.data.frame(do.call("rbind", cobn.list))
    names(all.orthologs) <- c("Efa", "Ete")
    all.orthologs$Ete <- gsub("Ete\\|", "", all.orthologs$Ete)
    all.orthologs$Efa <- gsub("Efa\\|", "", all.orthologs$Efa)
    return(all.orthologs)
}

all.orthologs <- get.all.pairedE.orthos()


get.all.pairedT.orthos <- function (){
    Efa.Tgo <- lapply(omcl, function (x) {
        x[grep("Efa\\||Tgo\\|", x)]
    })
    ## has to be a cluster still 
    Efa.Tgo <- Efa.Tgo[unlist(lapply(Efa.Tgo, length))>1]
    ## and contain genes from both
    Efa.Tgo <- Efa.Tgo[unlist(lapply(Efa.Tgo, function (x){ any(grepl("Efa\\|", x)) &
                                                                any(grepl("Tgo\\|", x)) }))]
    ## split them appart again
    Efa.list <-lapply(Efa.Tgo, function (x) x[grep("Efa\\|", x)])
    Tgo.list <- lapply(Efa.Tgo, function (x) x[grep("Tgo\\|", x)])
    ## Do all combinations of each Efa gene with with each Tgo gene in a cluster
    cobn.list <- lapply(1:length(Efa.list), function(i){
        exp <- rep(Efa.list[[i]], times=length(Tgo.list[[i]]))
        cbind(exp, Tgo.list[[i]])
    })
    all.orthologs <- as.data.frame(do.call("rbind", cobn.list))
    names(all.orthologs) <- c("Efa", "Tgo")
    all.orthologs$Tgo <- gsub("Tgo\\|", "", all.orthologs$Tgo)
    all.orthologs$Efa <- gsub("Efa\\|", "", all.orthologs$Efa)
    return(all.orthologs)
}

all.orthologs <- merge(all.orthologs, get.all.pairedT.orthos(), all.x=TRUE)

## read the Tenella data of Walker et al. 
Walker.data <- read.csv("data/12864_2015_1298_MOESM3_ESM.csv", skip=1)

names(Walker.data)[5:8] <- paste(names(Walker.data)[5:8], "GamVsMer", sep=".")
names(Walker.data) <- gsub("\\.1", "\\.GamVsSpor", names(Walker.data))

## make thins numeric
Walker.data <- Walker.data[, c("gene.ID", "gametocyte", "merozoite", "sporozoite",
                               "log2.fold.change.GamVsSpor", "log2.fold.change.GamVsMer")] 

Walker.data <- merge(all.orthologs, Walker.data, by.x="Ete", by.y="gene.ID")

Walker.RC <- Walker.data[, c("Efa", "Ete", "Tgo", "gametocyte", "merozoite", "sporozoite")]
names(Walker.RC) <-  c("Efa", "Ete", "Tgo", "Walker.Gam", "Walker.Mer", "Walker.Spo")

Walker.RC[,c("Walker.Gam", "Walker.Mer", "Walker.Spo")] <-
    log2(Walker.RC[,c("Walker.Gam", "Walker.Mer", "Walker.Spo")])
Walker.RC[is.na(Walker.RC)] <- 0

Walker.FC <-  Walker.data[, c("Efa", "Ete", "Tgo", "log2.fold.change.GamVsMer", "log2.fold.change.GamVsSpor")]

## More from Reid et al via ToxoDB (tedious manual downloads)
read.Reid.data <- function(){
    OOvsMer <- read.delim("data/OO_vs_Mer.txt")
    OOvsMer <- OOvsMer[, -ncol(OOvsMer)]
    names(OOvsMer) <- c("Ete", "FC.OOvsMer", "Reid.Mer", "Reid.Ooc")
    OOvsspOO <- read.delim("data/OO_vs_spOO.txt")
    OOvsspOO <- OOvsspOO[, -ncol(OOvsspOO)]
    names(OOvsspOO) <- c("Ete", "FC.OOvsspOO", "Reid.spOoc", "Reid.Ooc")
    OOvsSp <- read.delim("data/OO_vs_SP.txt")
    OOvsSp <- OOvsSp[, -ncol(OOvsSp)]
    names(OOvsSp) <- c("Ete", "FC.OOvsSp", "Reid.Sp", "Reid.Ooc")
    foo <- merge(OOvsMer, OOvsspOO, by="Ete")
    bar <- merge(foo, OOvsSp, by="Ete")
    bar <- bar[, !grepl("\\.x$|\\.y$", names(bar))]
    return(bar)
}

Reid.data <- read.Reid.data()
Reid.data <- merge(Reid.data, all.orthologs, by="Ete")

Reid.RC <- Reid.data[, c("Efa", "Ete", "Tgo", "Reid.Mer","Reid.spOoc", "Reid.Sp", "Reid.Ooc")]
Reid.FC <- Reid.data[, c("Efa", "Ete", "Tgo", "FC.OOvsMer","FC.OOvsspOO", "FC.OOvsSp")]

## Data from Hehl lab on toxo in the cat
read.Hehl.data <- function(){
    Tg7vs3 <- read.delim("data/Tg7vs3.txt")
    Tg7vs3 <- Tg7vs3[, -ncol(Tg7vs3)]
    names(Tg7vs3) <- c("Tgo", "FC.Tg7vs3", "Hehl.Tg3", "Hehl.Tg7", "prev")
    Tg7vsTach <- read.delim("data/Tg7vsTach.txt")
    ## remove also the previous id column and the unnecessary product col
    Tg7vsTach <- Tg7vsTach[, -c((ncol(Tg7vsTach)-2):ncol(Tg7vsTach))]
    names(Tg7vsTach) <- c("Tgo", "FC.Tg7vsTach", "Hehl.TgTach", "Hehl.Tg7")
    bar <- merge(Tg7vs3, Tg7vsTach, by="Tgo")
    bar <- bar[, !grepl("\\.y$", names(bar))]
    foo <- strsplit(as.character(bar$prev), ", ")
    foobar <- lapply(foo, function (x) grep("^TGME49_", x, value=TRUE))
    ### discard where we have multiple or no previous ids mathich our ids
    bar <- bar[unlist(lapply(foobar, function (x) length(x)==1)), ]
    foobar <- foobar[unlist(lapply(foobar, function (x) length(x)==1))]
    bar$Tgo <- unlist(foobar)
    return(bar[, !names(bar)%in%"prev"])
}

Hehl.data <- read.Hehl.data()
Hehl.data <- merge(all.orthologs, Hehl.data, by="Tgo")

Hehl.RC <- Hehl.data[, c("Efa", "Ete", "Tgo", "Hehl.Tg3", "Hehl.Tg7.x", "Hehl.TgTach")]
Hehl.FC <- Hehl.data[, c("Efa", "Ete", "Tgo", "FC.Tg7vs3", "FC.Tg7vsTach")]

## Work with normalized counts
get.ortholog.RC <- function(){
    Ef.counts <- log2(cpm(Ef.1st.pass.model[[4]]))
    Ef.counts[is.infinite(Ef.counts)] <- 0
    Ef.counts <- as.data.frame(Ef.counts)
    Ef.counts$Efa <- rownames(Ef.counts)
    foo <- merge(Hehl.RC, Ef.counts)
    bar <- merge(Reid.RC, Walker.RC)
    bar <- unique(bar)
    foobar <- merge(bar, foo)
    foobar[, 3:ncol(foobar)] <- apply(foobar[, 3:ncol(foobar)], 2,
                                      function (x) as.numeric(as.character(x)))
    return(foobar)
}

RC.ortho <- merge(RC.ortho, cpm(Ef.1st.pass.model[[4]]),
                  by.x="Efa", by.y=0)
RC.ortho[,3:ncol(RC.ortho)] <- apply(RC.ortho[,3:ncol(RC.ortho)], 2,
                                     function (x) as.numeric(as.character(x)))
=======
RC.ortho <- get.ortholog.RC()

CS <- colSums(RC.ortho[,4:ncol(RC.ortho)])
scaling.f <- max(CS)/CS
RC.ortho[, 4:ncol(RC.ortho)] <- t(t(RC.ortho[, 4:ncol(RC.ortho)])*scaling.f)

pdf("figuresANDmanuscript/RC_correlation_Efal_Ete_Tgo.pdf", heigh=15, width=15)
pheatmap(cor(RC.ortho[,4:ncol(RC.ortho)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()

###
mean.columns <- function(x){
  reps <- as.factor(gsub("_rep\\d", "", names(x)))
  y <- do.call(rbind, by(t(x), reps, colMeans))
  t(y)
}
  
RC.ortho.mean <- as.data.frame(cbind(RC.ortho[, 1:13],
                                     mean.columns(RC.ortho[, 14:ncol(RC.ortho)])))

pdf("figuresANDmanuscript/RC_mean_correlation_Efal_Ete_Tgo.pdf", width=10, height=10)
pheatmap(cor(RC.ortho.mean[,4:ncol(RC.ortho.mean)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()


png("figuresANDmanuscript/RC_mean_correlation_Efal_Ete_Tgo.png", width=1500, height=1500, res=150)
pheatmap(cor(RC.ortho.mean[,4:ncol(RC.ortho.mean)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()



non.one.to.one <- RC.ortho.mean$Efa[duplicated(RC.ortho.mean$Efa)]
RC.ortho.mean.oneONE <- RC.ortho.mean[!RC.ortho$Efa%in%non.one.to.one, ]

pdf("figuresANDmanuscript/RC_mean_correlation_Efal_Ete_Tgo_ONE_ONE.pdf", width=10, height=10)
pheatmap(cor(RC.ortho.mean.oneONE[,4:ncol(RC.ortho.mean.oneONE)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()


## limit this to lifecycle significant genes?!
union.of.Ef.cycle.diff.top.100 <-
    unique(unlist(lapply(Ef.1st.pass.model[[3]][1:10], head, n=500)))

RC.ortho.mean.oneONE.LCsig <-
    RC.ortho.mean.oneONE[RC.ortho.mean.oneONE$Efa%in%union.of.Ef.cycle.diff.top.100, ]

pdf("figuresANDmanuscript/RC_mean_correlation_Efal_Ete_Tgo_ONE_ONE_LCsig.pdf", width=10, height=10)
pheatmap(cor(RC.ortho.mean.oneONE.LCsig[,4:ncol(RC.ortho.mean.oneONE.LCsig)],
             method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()

## ### Work with Fold Changes... ## What to do wiht Max and Min? ## probably no good idea
## FC.ortho <- merge(all.orthologs,
##                   TenGam[, grepl("gene.ID|log2.fold.change", names(TenGam))],
##                   by.x="Ete", by.y="gene.ID")
## FC.ortho <- FC.ortho[!rowSums(is.na(FC.ortho))>0, ]

## all.Ef.FC.list <- lapply(Ef.1st.pass.model[[2]], function(x) {
##   as.data.frame(cbind(gene=rownames (x), logFC=x[,1]))
## })

## all.FC.df <- Reduce(function (x,y)merge(x,y, by="gene"),all.Ef.FC.list)
## names(all.FC.df) <- c("gene.ID", paste("logFC", names(Ef.1st.pass.model[[2]]), sep="."))

## FC.ortho <- merge(FC.ortho,
##                   all.FC.df,
##                   by.x="Efa", by.y="gene.ID")

## FC.ortho[,3:ncol(FC.ortho)] <- apply(FC.ortho[,3:ncol(FC.ortho)],
##                                      2, function(x) as.numeric(as.character(x)))

## pdf("figuresANDmanuscript/TenellaGamVsMer_vs_N3vsN7.pdf")
## ggplot(FC.ortho, aes(log2.fold.change.GamVsMer, logFC.N3vsN7*-1)) +
##   geom_point(alpha=0.8) +
##   stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
##   scale_fill_gradient(low = "yellow", high = "red") +
##   scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
##   stat_smooth() +
##   theme_bw()
## dev.off()


## pdf("figuresANDmanuscript/TenellaGamVsMer_vs_N5vsN7.pdf")
## ggplot(FC.ortho, aes(log2.fold.change.GamVsMer, logFC.N5vsN7*-1)) +
##   geom_point(alpha=0.8) +
##   stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
##   scale_fill_gradient(low = "yellow", high = "red") +
##   scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
##   stat_smooth() +
##   theme_bw()
## dev.off()


## pdf("figuresANDmanuscript/TenellaGamVsSpor_vs_N7vsSpor.pdf")
## ggplot(FC.ortho, aes(log2.fold.change.GamVsSpor, logFC.N7vsSp)) +
##   geom_point(alpha=0.8) +
##   stat_density2d(aes(alpha=..level.., fill=..level..), size=2, geom="polygon") +
##   scale_fill_gradient(low = "yellow", high = "red") +
##   scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
##   stat_smooth() +
##   theme_bw()
## dev.off()

