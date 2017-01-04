library(ggplot2)
library(pheatmap)
library(reshape)
library(ggplot2)
library(gridExtra)
library(RSvgDevice)
library(xtable)

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

  #Haemosporidia.clusters <- rownames(ortholog.presence[select.deeper.nodes(Haemosporidia, ), ])
  ## with or without comma?????????
  
  Haemosporidia.clusters <- rownames(ortholog.presence[select.deeper.nodes(Haemosporidia), ])
  Coccidia.clusters <- rownames(ortholog.presence[select.deeper.nodes(Coccidia), ])
  Coccidia.clusters <- Coccidia.clusters[!Coccidia.clusters%in%Eimeria.clusters &
                                         !Coccidia.clusters%in%Sarcocystida.clusters]

  Coccidia.genes <- get.genes.4.clusters(omcl[Coccidia.clusters])

  Aconoidasida.clusters <- rownames(ortholog.presence[select.deeper.nodes(Aconoidasida), ])
  Aconoidasida.clusters <- Aconoidasida.clusters[!Aconoidasida.clusters%in%Piroplasma.clusters &
                                                 !Aconoidasida.clusters%in%Haemosporidia.clusters]

  Apicomplexa.clusters <- rownames(ortholog.presence[select.deeper.nodes(Apicomplexa), ])
  Apicomplexa.clusters <- Apicomplexa.clusters[!Apicomplexa.clusters%in%Aconoidasida.clusters &
                                               !Apicomplexa.clusters%in%Coccidia.clusters &
                                               !Apicomplexa.clusters%in%Eimeria.clusters &
                                               !Apicomplexa.clusters%in%Ef.only.clusters]
  Apicomplexa.genes <- get.genes.4.clusters(omcl[Apicomplexa.clusters])

  ApicomplexaC.clusters <- rownames(ortholog.presence[select.deeper.nodes(ApicomplexaC), ])
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
    Haemosporidi=Haemosporidia.clusters, # 
    Piroplasma=Piroplasma.clusters,     #  Aplicomplexan parasite
    Aconoidasida=Aconoidasida.clusters, # Aplicomplexan parasite
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
hcluster <- list()
hcluster[["Mm"]] <- read.table("output_data/Mm_hclustered_cycle.csv", sep=",",
                               header=TRUE)
hcluster[["Ef"]] <- read.table("output_data/Ef_hclustered_cycle.csv", sep=",",
                               header=TRUE)
###


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
lots.of.f$FDR <- p.adjust(lots.of.f$value.p.value, method="BH")

lots.signif <- lots.of.f[ lots.of.f[, "FDR"]<0.01, ]
names(lots.signif) <- c("Ef_Cluster", "conservation", "p-value", "odds.ratio", "FDR")

tabl.lots <-
    xtable(lots.signif[, c("Ef_Cluster", "conservation",
                           "odds.ratio", "p-value", "FDR")],
           digits=c(NA, NA, NA, 2, -2, -2))

print(tabl.lots,
      type = "html", file = "tables/Table3_Conserv_Cluster.html", include.rownames = F,
      format.args = list(big.mark = ",", decimal.mark = "."))



#################################################
## Ortholog expression analysis############
#########################################

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

## downloaded from ToxoDB v. 28
## Walker data
Walker.RC <- read.delim("data/Walker_profiles.min.pct")
names(Walker.RC) <- c("Ete", "Walker.Gametocytes",
                      "Walker.Merozoites", "Walker.Sporozoites")
Walker.RC <- merge(Walker.RC, all.orthologs)


## More from Reid et al via ToxoDB 28
Reid.RC <- read.delim("data/Reid_profiles.min.pct")
names(Reid.RC)[1] <- "Ete"
names(Reid.RC)[2:ncol(Reid.RC)] <-
    paste("Reid", names(Reid.RC)[2:ncol(Reid.RC)], sep=".")
Reid.RC <- merge(Reid.RC, all.orthologs)


Hehl.RC <- read.delim("data/Hehl_profiles.min.pct")
names(Hehl.RC)[1] <- "Tgo_Now"
names(Hehl.RC)[2:ncol(Hehl.RC)] <-
    paste("Hehl", names(Hehl.RC)[2:ncol(Hehl.RC)], sep=".")

ID.toxo <- read.delim("data/ID_mapping_toxo.txt")
Tg.prev <- strsplit(as.character(ID.toxo$X.Previous.ID.s..), ", ")
Tg.prev <- lapply(Tg.prev, function (x) grep("^TGME49_\\d", x, value=TRUE))
Tg.prev <- lapply(Tg.prev, "[", 1)
Tg.prev <- unlist(Tg.prev)
Tgo <- gsub("<br>", "", Tg.prev)
Tg.link <- as.data.frame(cbind(as.character(ID.toxo$X.Gene.ID.), Tgo))

Hehl.RC <- merge(Hehl.RC, Tg.link, by.x = "Tgo_Now", by.y = "V1")
Hehl.RC$Tgo_Now <- NULL
Hehl.RC <- Hehl.RC[!is.na(Hehl.RC$Tgo),]

Hehl.RC <- merge(Hehl.RC, all.orthologs)

## Work with normalized counts
get.ortholog.RC <- function(){
    Ef.counts <- log2(Ef.nRC)
    Ef.counts[is.infinite(Ef.counts)] <- 0
    Ef.counts <- as.data.frame(Ef.counts)
    Ef.counts$Efa <- rownames(Ef.counts)
    foo <- merge(Hehl.RC, Ef.counts)
    foo <- merge(foo, Walker.RC)
    foo <- merge(foo, Reid.RC)
    bar <- unique(foo)
    bar[, 3:ncol(bar)] <- apply(bar[, 3:ncol(bar)], 2,
                                function (x) as.numeric(as.character(x)))
    return(bar)
}

RC.ortho <- get.ortholog.RC()


mean.columns <- function(x){
  reps <- as.factor(gsub("_rep\\d", "", names(x)))
  y <- do.call(rbind, by(t(x), reps, colMeans))
  t(y)
}

RC.ortho.mean <- as.data.frame(cbind(RC.ortho[, 1:13],
                                     mean.columns(RC.ortho[, 14:ncol(RC.ortho)])))


non.one.to.one <- RC.ortho.mean$Efa[duplicated(RC.ortho.mean$Efa)]

RC.ortho.oneONE <- RC.ortho[!RC.ortho$Efa%in%non.one.to.one, ]

RC.ortho.mean.oneONE <- RC.ortho.mean[!RC.ortho.mean$Efa%in%non.one.to.one, ]


pdf("figures/Figure4_RC_mean_correlation_Efal_Ete_Tgo_ONEONE.pdf", width=10, height=10, onefile=FALSE)

## The figure
#devSVG("figures/Figure4a_RC_mean_correlation_Efal_Ete_Tgo_ONEONE.svg", width=10, height=10, onefile=FALSE)
pheatmap(cor(RC.ortho.mean.oneONE[,4:ncol(RC.ortho.mean.oneONE)], method="spearman"),
         display_numbers=TRUE, scale="none")
dev.off()


