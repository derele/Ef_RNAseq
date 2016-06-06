library(ggplot2)
library(reshape)
library(plyr)
library(phangorn)
#library(multicore)

## Orthomcl
omcl <- readLines("data/Omcl_groups.txt")
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

all.orthologs <- do.call("rbind", cobn.list)
