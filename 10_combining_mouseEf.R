library(ggplot2)
library(pheatmap)
library(doParallel) # more flexible than the library parallel
library(RSvgDevice)

## raw count data
Mm.nRC <- read.table("output_data/Mm_norm_counts.csv", sep=",")
Ef.nRC <- read.table("output_data/Ef_norm_counts.csv", sep=",")

both.col <- intersect(colnames(Mm.nRC), colnames(Ef.nRC))

## Reid and Berriman approach: We randomized all profiles and
## calculated the Pearson correlation coefficient between every
## interspecific pair of genes. We repeated this 10^5 times and
## calculated the P-value as the number of times a randomized pair of
## genes profiles were correlated at least as well as the real
## profiles, divided by 10^5

## Cors <- cor(t(Mm.nRC[, both.col]), t(Ef.nRC[, both.col]))

## registerDoParallel(cores=20)

## outer.reps <- 5000
## inner.reps <- 100

## ## sum up parallel the numbers the random correlation values are
## ## bigger than the real correlations
## parallel.RnB <- function(inner.reps){
##     randCors <-
##         foreach(icount(inner.reps), .combine="+") %dopar% {
##             rand.col <- sample(both.col, length(both.col))
##             RC <- cor(t(Mm.nRC[, rand.col]),
##                       t(Ef.nRC[, both.col]))        
##             res <- abs(RC)>abs(Cors)
##             ## avoid an accumulation of NAs for In cor(t(Mm.nRC[,
##             ## both.col]), t(Ef.nRC[, both.col])) : the standard
##             ## deviation is zero
##             res[is.na(res)] <- TRUE
##             return(res)
##         }
##     return(randCors)
## }


## ## a large number of non-parallelized chunks. Putting everything in
## ## the parallel function would eat up memory
## RnB <- matrix(0, nrow=nrow(Cors), ncol=ncol(Cors))
## for(i in 1:outer.reps){
##     ## sum each result with what has been there before
##     RnB.new <- parallel.RnB(inner.reps)
##     print(paste("finished outer rep", i))
##     RnB <- RnB + RnB.new
## }

## file.RnB <- paste("//home/ele/RnB_",
##                   round(as.vector(Sys.time())), ".Rdata", sep="")

## save(RnB, file = file.RnB)


## RnB.final <- RnB / (outer.reps*inner.reps)

## prod.file.RnB <- paste("/home/ele/RnB_Prod_",
##                        round(as.vector(Sys.time())), ".Rdata", sep="")

## save(RnB.final, file = prod.file.RnB)

load("/SAN/Eimeria_Totta/RnB_1478263654.Rdata")
load("/SAN/Eimeria_Totta/RnB_Prod_1478263755.Rdata")

is.zero <- (RnB.final==0)

is.zero.all.cols <- apply(is.zero, 2, any)
Ef.genes.ineracting <- colnames(RnB)[is.zero.all.cols] ## TYPO interacting/ineracting

is.zero.all.rows <- apply(is.zero, 1, any)
Mm.genes.ineracting <- rownames(RnB)[is.zero.all.rows]

## the interaction network for later
inter.net <- is.zero[is.zero.all.rows, is.zero.all.cols]

interactions.per.Ef.gene <- data.frame(apply(is.zero, 2, sum))
names(interactions.per.Ef.gene) <- "num.interactions"
interactions.per.Ef.gene$from.to <- "Ef.to.Mm"

interactions.per.Mm.gene <- data.frame(apply(is.zero, 1, sum))
names(interactions.per.Mm.gene) <- "num.interactions"
interactions.per.Mm.gene$from.to <- "Mm.to.Ef"

interactions.per.gene <- rbind(interactions.per.Ef.gene,
                               interactions.per.Mm.gene)


## plot mouse and parasite densities together
devSVG(file = "figures/Figure5a_dis.svg", height = 5, width = 6)
ggplot(interactions.per.gene, aes(x=num.interactions, y=..count..,
                                  color=from.to))+
    geom_freqpoly(stat="bin", binwidth=1)+
    scale_y_sqrt("Times reported", breaks=seq(0, 100, by=10)^2)+
    scale_x_sqrt("Number of interactions", breaks=seq(0, 35, by=5)^2) +
    ggtitle("Interactions per gene") +
    theme_bw()
dev.off()


library(ggnet)
library(sna)
library(intergraph)


net <- network(inter.net, directed = FALSE,
               matrix.type = "bipartite")

col <- c("actor" = "red", "event" = "blue")

pdf("figures/test.pdf")
ggnet2(net, color = "mode", palette = col, size=0.5,
       edge.alpha = 0.5, edge.width = 0.5)
dev.off()


rownames(hcluster[["Ef"]])

pdf("figures/test.pdf")
ggnet2(net, color = "mode", palette = col, size=0.5,
       edge.alpha = 0.5, edge.width = 0.5)
dev.off()
