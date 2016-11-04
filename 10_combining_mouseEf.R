library(ggplot2)
library(pheatmap)
library(doParallel) # more flexible than the library parallel

## raw count data
Mm.nRC <- read.table("output_data/Mm_norm_counts.csv", sep=",")
Ef.nRC <- read.table("output_data/Ef_norm_counts.csv", sep=",")

both.col <- intersect(colnames(Mm.nRC), colnames(Ef.nRC))

## nRC <- rbind(Mm.nRC[, both.col], Ef.nRC[, both.col])

## get.scaled.and.clustered <- function(data){
##     rMeans <- rowMeans(data)
##     rSD <- apply(data, 1, sd)
##     scaled.data <- scale(t(data), rMeans, rSD)
##     ## NAs for whenever tehere is no variance/sd need to be removed
##     scaled.data <- scaled.data [, apply(scaled.data, 2, function(x) all(!is.na(x) ))]
##     hclustered <- hclust(as.dist(t(1-abs(cor(scaled.data, method="spearman")))))
##     return(hclustered)
## }

## nRC.clusters <- get.scaled.and.clustered(nRC)

## nRC.cuttree <- as.data.frame(cutree(nRC.clusters, k=1000))
## names(nRC.cuttree) <- "cluster"
## nRC.cuttree$cluster <- paste("CL", nRC.cuttree$cluster, sep="_")

## nRC.cuttree$species <- ifelse(grepl("^ENSMUS", rownames(nRC.cuttree)), "Mm", "Ef")

## spec.no <- tapply(nRC.cuttree$species, nRC.cuttree$cluster, function(x) {
##     cbind(Mm=length(x[x%in%"Mm"]),
##           Ef=length(x[x%in%"Ef"]))
## })


## sp.per.cluster <- as.data.frame(do.call(rbind, spec.no))
## sp.per.cluster$cluster <- paste("CL", 1:nrow(boo), sep="_")

## ggplot(sp.per.cluster, aes(Mm, Ef)) + geom_jitter() #+ scale_x_log10() + scale_y_log10()

## ## arbitrary thresholds for coexpression
## sp.per.cluster.co <- sp.per.cluster[sp.per.cluster$Mm>5&sp.per.cluster$Ef>3,]

## nRC.cuttree.co <- nRC.cuttree[nRC.cuttree$cluster%in%sp.per.cluster.co$cluster, ]

## nRC.co <- nRC[rownames(nRC)%in%rownames(nRC.cuttree.co), ]

## pheatmap(nRC.co, scale="row",
##          annotation_row = nRC.cuttree.co,
##          cutree_rows = 7, 
##          show_rownames = F
##          )

## Reid and Berriman approach: We randomized all profiles and
## calculated the Pearson correlation coefficient between every
## interspecific pair of genes. We repeated this 10^5 times and
## calculated the P-value as the number of times a randomized pair of
## genes profiles were correlated at least as well as the real
## profiles, divided by 10^5

Cors <- cor(t(Mm.nRC[, both.col]), t(Ef.nRC[, both.col]))

registerDoParallel(cores=20)

outer.reps <- 5000
inner.reps <- 100

## uncomment for testing
## outer.reps <- 10
## inner.reps <- 20

## sum up parallel the numbers the random correlation values are
## bigger than the real correlations
parallel.RnB <- function(inner.reps){
    randCors <-
        foreach(icount(inner.reps), .combine="+") %dopar% {
            rand.col <- sample(both.col, length(both.col))
            RC <- cor(t(Mm.nRC[, rand.col]),
                      t(Ef.nRC[, both.col]))        
            res <- abs(RC)>abs(Cors)
            ## avoid an accumulation of NAs for In cor(t(Mm.nRC[,
            ## both.col]), t(Ef.nRC[, both.col])) : the standard
            ## deviation is zero
            res[is.na(res)] <- TRUE
            return(res)
        }
    return(randCors)
}


## a large number of non-parallelized chunks. Putting everything in
## the parallel function would eat up memory
RnB <- matrix(0, nrow=nrow(Cors), ncol=ncol(Cors))
for(i in 1:outer.reps){
    ## sum each result with what has been there before
    RnB.new <- parallel.RnB(inner.reps)
    print(paste("finished outer rep", i))
    RnB <- RnB + RnB.new
}


file.RnB <- paste("//home/ele/RnB_",
                  round(as.vector(Sys.time())), ".Rdata", sep="")

save(RnB, file = file.RnB)


RnB.final <- RnB / (outer.reps*inner.reps)

prod.file.RnB <- paste("/home/ele/RnB_Prod_",
                       round(as.vector(Sys.time())), ".Rdata", sep="")

save(RnB.final, file = prod.file.RnB)

## all.x <- apply(RnB.final, 2, function(x) sum(x==0))

foo <- (RnB.final==0)

bar <- apply(foo, 2, any)
Ef.genes.ineracting <- colnames(RnB)[bar]

baz <- apply(foo, 1, any)
Mm.genes.ineracting <- rownames(RnB)[baz]


