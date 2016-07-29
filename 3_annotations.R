## needed from previous script 2_edgeR_diff.R:
## a TopTags (edgeR) object listing genes tested for differential
## expression

Mm.DE.test <- read.table("output_data/Mm_DEtest.csv", sep=",")
Ef.DE.test <- read.table("output_data/Ef_DEtest.csv", sep=",")


## The Bioconductor mouse annotation libary
library(Mus.musculus)

## merge all this information in a data frame to combine mouse
## annotations with gene_id from Cufflinks get the information from
## this annotation dbi into data frames
.get.annot.frame <- function(){
    ens2geneID <- as.data.frame(org.Mm.egENSEMBL)
    name2SY <- as.data.frame(org.Mm.egSYMBOL2EG) # get 5 digit "symbol" to "short gene name" mapping
    SY2Lname<- as.data.frame(org.Mm.egGENENAME) # get 5 digit "symbol" to "long name name" mapping
    SY2GO <- as.data.frame(org.Mm.egGO) # get 5 digit "symbol" to gene ontology mapping
    annot.frame <- merge(name2SY, ens2geneID, by = "gene_id")
    annot.frame <- merge(annot.frame, SY2GO, by = "gene_id")
    annot.frame <- merge(annot.frame, SY2Lname, by = "gene_id")
    return(annot.frame)
}

annot.frame <- .get.annot.frame()

## ensembl2GO contains the mapping of ensembllink ids to GO terms
gene2GO <- list()

gene2GO[["Mm"]] <- by(annot.frame, annot.frame$ensembl_id,
                      function(x) as.character(x$go_id))
## define the universe of tested genes
gene2GO[["Mm"]] <- gene2GO[["Mm"]][names(gene2GO[["Mm"]])%in%
                                   unique(Mm.DE.test$gene)]


go_Ef <- read.csv("data/12864_2014_6777_GO-annot.csv")

gene2GO[["Ef"]] <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))
gene2GO[["Ef"]] <- gene2GO[["Ef"]][names(gene2GO[["Ef"]])%in%
                                   unique(Ef.DE.test$gene)]

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

#####  objects created here for further use:
## 1. annot.frame - a data frame containing coprehensive annotation information
## 2. gene2GO - a list linking cufflinks ids to GO terms for both mouse and Eimeria

## 3. gene.cluster.category - a categorization of genes, gene families
## (clusters) into phylogenetic categories of "novelty"

