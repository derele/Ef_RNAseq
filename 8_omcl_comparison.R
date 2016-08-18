library(ggplot2)
library(pheatmap)
library(reshape)
library(ggplot2)

Ef.nRC <- read.table("output_data/Ef_norm_counts.csv", sep=",", header=TRUE)
Ef.nRC <- matrix(unlist(Ef.nRC), nrow=nrow(Ef.nRC), ncol=ncol(Ef.nRC),
                 dimnames=list(rownames(Ef.nRC), colnames(Ef.nRC)))

### ORTHOMCL categories
read.omcl <- function(path="data/Omcl_groups.txt"){
  omcl <- readLines(path)
  omcl <- lapply(omcl, strsplit, " ")
  omcl <- lapply(omcl, function(x) x[[1]])
  names(omcl) <- lapply(omcl, function(x) gsub("\\:", "", x[1]))
  omcl <- lapply(omcl, function(x) x[-1])

  ## make gene-ID identifiers used in the omcl clsters
  omcl <- lapply(omcl, function (x)
                 gsub("(Efa|)NODE_(.*_\\d+.g\\d+).t\\d+","\\1EfaB_\\2",  x))

  omcl <- lapply(omcl, function (x) x[!duplicated(x)])
  omcl <- omcl[!unlist(lapply(omcl, length))==1]
  omcl <- omcl[!names(omcl)%in%"api10278"] ## exclude a cluster... wonder now why that ;-)

}

omcl <- read.omcl()

get.genes.4.clusters <- function (clusters) {
  gsub("Efa\\|", "",   grep("Efa\\|",
                            unlist(clusters),
                            value=TRUE))
}


get.omcl.categories <- function(ortho.obj){
  ## Orthomcl
  get.presence.omcl <- function(ortho.obj){
    omcl.tab <- lapply(ortho.obj, function (x) table(substr(x, 1, 3)))
    spp.inf <- unique(unlist(lapply(omcl.tab, names)))
    dat <- data.frame()
    for(i in seq(along=omcl.tab)) for(j in names(omcl.tab[[i]])){
      dat[i,j] <- omcl.tab[[i]][j]}
    rownames(dat) <- names(omcl.tab)
    dat[is.na(dat)] <- 0
    ortholog.presence <- as.data.frame(apply(dat, 2, function (x) as.numeric(x>0)))
    rownames(ortholog.presence) <- rownames(dat)
    return(ortholog.presence)
  }

  ortholog.presence <- get.presence.omcl(omcl)

  Ef.only.clusters <- names(omcl[ortholog.presence$Efa==1 & rowSums(ortholog.presence)==1 ])
  Ef.only.genes <- get.genes.4.clusters(omcl[Ef.only.clusters])

  Eimeria <- c("Efa", "Ema", "Ete")
  Sarcocystida <- c("Tgo", "Nca")
  Coccidia <- c(Eimeria, Sarcocystida)
  Piroplasma <- c("Bbo", "Tan")
  Haemosporidia <- c("Pfa", "Pbe")
  Aconoidasida <- c(Piroplasma, Haemosporidia)
  Apicomplexa <- c(Coccidia, Aconoidasida)
  ApicomplexaC <- c(Apicomplexa, "Cho")

  ## use x=1 to restrict to multiple species clusters at terminal nodes 
  select.deeper.nodes <- function (Node, x=1){
    rowSums(ortholog.presence[, Node])>x & 
      rowSums(ortholog.presence[,!names(ortholog.presence)
                                %in%Node])==0
  }

  Eimeria.clusters <- rownames(ortholog.presence[select.deeper.nodes(Eimeria), ])
  Eimeria.clusters <- Eimeria.clusters[!Eimeria.clusters%in%Ef.only.clusters]

  Eimeria.genes <- get.genes.4.clusters(omcl[Eimeria.clusters])

  Sarcocystida.clusters <- rownames(ortholog.presence[select.deeper.nodes(Sarcocystida), ])

  Piroplasma.clusters <- rownames(ortholog.presence[select.deeper.nodes(Piroplasma), ])

  Haemosporidia.clusters <- rownames(ortholog.presence[select.deeper.nodes(Haemosporidia, ), ])

  Coccidia.clusters <- rownames(ortholog.presence[select.deeper.nodes(Coccidia, ), ])
  Coccidia.clusters <- Coccidia.clusters[!Coccidia.clusters%in%Eimeria.clusters &
                                         !Coccidia.clusters%in%Sarcocystida.clusters]

  Coccidia.genes <- get.genes.4.clusters(omcl[Coccidia.clusters])

  Aconoidasida.clusters <- rownames(ortholog.presence[select.deeper.nodes(Aconoidasida, ), ])
  Aconoidasida.clusters <- Aconoidasida.clusters[!Aconoidasida.clusters%in%Piroplasma.clusters &
                                                 !Aconoidasida.clusters%in%Haemosporidia.clusters]

  Apicomplexa.clusters <- rownames(ortholog.presence[select.deeper.nodes(Apicomplexa, ), ])
  Apicomplexa.clusters <- Apicomplexa.clusters[!Apicomplexa.clusters%in%Aconoidasida.clusters &
                                               !Apicomplexa.clusters%in%Coccidia.clusters &
                                               !Apicomplexa.clusters%in%Eimeria.clusters &
                                               !Apicomplexa.clusters%in%Ef.only.clusters]
  Apicomplexa.genes <- get.genes.4.clusters(omcl[Apicomplexa.clusters])

  ApicomplexaC.clusters <- rownames(ortholog.presence[select.deeper.nodes(ApicomplexaC,), ])
  ApicomplexaC.clusters <- ApicomplexaC.clusters[!ApicomplexaC.clusters%in%Apicomplexa.clusters&
                                                 !ApicomplexaC.clusters%in%Coccidia.clusters &
                                                 !ApicomplexaC.clusters%in%Eimeria.clusters &
                                                 !ApicomplexaC.clusters%in%Ef.only.clusters]

  ApicomplexaC.genes <- get.genes.4.clusters(omcl[ApicomplexaC.clusters])


  cluster.grouping <- melt(list(
    Efalciformis=Ef.only.clusters,
    Eimeria=Eimeria.clusters,
    Sarcocystida=Sarcocystida.clusters,
    Coccidia=Coccidia.clusters,
    Haemosporidi=Haemosporidia.clusters,
    Piroplasma=Piroplasma.clusters,
    Aconoidasida=Aconoidasida.clusters,
    Apicomplexa=Apicomplexa.clusters,
    ApicomplexaC=ApicomplexaC.clusters))

  names(cluster.grouping) <- c("cluster", "grouping")
  
  Ef.gene.ortho.list <- lapply(omcl, get.genes.4.clusters)
  Ef.gene.ortho.list <- Ef.gene.ortho.list[unlist(lapply(Ef.gene.ortho.list, function (x) length(x)>0))]
  Ef.gene.ortho.df <- melt(Ef.gene.ortho.list)
  names(Ef.gene.ortho.df) <- c("gene", "cluster")

  gene.clus.cat <- merge(Ef.gene.ortho.df, cluster.grouping, all.x=TRUE)
  gene.clus.cat$grouping[is.na(gene.clus.cat$grouping)] <- "conserved"

  return(gene.clus.cat)
}

gene.cluster.category <- get.omcl.categories(omcl)

write.table(gene.cluster.category, "output_data/gene_family_category.csv", sep=",")


## testing for overrepresentation in conservating grouping


lots.of.f <- lapply(unique(gene.cluster.category$grouping), function (cons.group){
    inner <-
        lapply(unique(hcluster[["Ef"]]$Cluster), function (hier.cluster){
            fisher.test(rownames(Ef.nRC)%in%
                        gene.cluster.category[gene.cluster.category$grouping%in%cons.group,"gene"], 
                        rownames(Ef.nRC)%in%
                        rownames(hcluster[["Ef"]])[hcluster[["Ef"]]$Cluster%in%hier.cluster])[c("p.value",
                                                                                                "estimate")]
        })
    names(inner) <- paste("Cluster", unique(hcluster[["Ef"]]$Cluster), sep="_")
    return(inner)
})

names(lots.of.f) <- unique(gene.cluster.category$grouping)

lots.of.f <- melt(lots.of.f)
lots.of.f <- reshape(lots.of.f, idvar = c("L2", "L1"), timevar="L3", direction="wide")
lots.of.f$FDR <- round(p.adjust(lots.of.f$value.p.value, method="BH"), 4)

lots.of.f[ lots.of.f[, "FDR"]<0.01, ]

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

Walker.data <- Walker.data[, c("gene.ID", "gametocyte", "merozoite", "sporozoite")]

Walker.data <- merge(all.orthologs, Walker.data, by.x="Ete", by.y="gene.ID")

Walker.RC <- Walker.data[, c("Efa", "Ete", "Tgo", "gametocyte", "merozoite", "sporozoite")]
names(Walker.RC) <-  c("Efa", "Ete", "Tgo", "Walker.Gam", "Walker.Mer", "Walker.Spo")

## log2 transform to make it compatible with ToxoDB
Walker.RC[,c("Walker.Gam", "Walker.Mer", "Walker.Spo")] <-
    log2(Walker.RC[,c("Walker.Gam", "Walker.Mer", "Walker.Spo")])

Walker.RC[is.na(Walker.RC)] <- 0

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

## Work with normalized counts
get.ortholog.RC <- function(){
    Ef.counts <- log2(Ef.nRC)
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

RC.ortho <- get.ortholog.RC()
RC.ortho[is.na(RC.ortho)] <- 0

pdf("Supplement/RC_correlation_Efal_Ete_Tgo_MESSED.pdf", heigh=15, width=15, onefile=FALSE)
pheatmap(cor(RC.ortho[,4:ncol(RC.ortho)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()

non.one.to.one <- RC.ortho.mean$Efa[duplicated(RC.ortho.mean$Efa)]

RC.ortho.oneONE <- RC.ortho[!RC.ortho$Efa%in%non.one.to.one, ]

pdf("Supplement/RC_correlation_Efal_Ete_Tgo_ONEONE.pdf", heigh=15, width=15, onefile=FALSE)
pheatmap(cor(RC.ortho.oneONE[,4:ncol(RC.ortho.oneONE)], method="spearman"),
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

RC.ortho.mean.oneONE <- RC.ortho.mean[!RC.ortho.mean$Efa%in%non.one.to.one, ]

pdf("Supplement/RC_mean_correlation_Efal_Ete_Tgo_MESSED.pdf", width=10, height=10, onefile=FALSE)
pheatmap(cor(RC.ortho.mean[,4:ncol(RC.ortho.mean)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()

pdf("figures/Figure4a_RC_mean_correlation_Efal_Ete_Tgo_ONEONE.pdf", width=10, height=10, onefile=FALSE)
pheatmap(cor(RC.ortho.mean.oneONE[,4:ncol(RC.ortho.mean.oneONE)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()


