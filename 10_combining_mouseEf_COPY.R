library(ggplot2)
library(pheatmap)
library(doParallel) # more flexible than the library parallel

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

## nRC <- rbind(Mm.nRC[, both.col], Ef.nRC[, both.col])

## get.scaled.and.clustered <- function(data){
##     rMeans <- rowMeans(data)
##     rSD <- apply(data, 1, sd)
##     scaled.data <- scale(t(data), rMeans, rSD)
##     ## NAs for whenever there is no variance/sd need to be removed
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

######## producing visual output
## make df:s with Ef and Mm data
interactions.per.Ef.gene <- data.frame(apply(is.zero, 2, sum))
interactions.per.Mm.gene <- data.frame(apply(is.zero, 1, sum))

interactions.per.Ef.gene <- tibble::rownames_to_column(interactions.per.Ef.gene, "gene")
interactions.per.Mm.gene <- tibble::rownames_to_column(interactions.per.Mm.gene, "gene")

colnames(interactions.per.Ef.gene)[2] <- "Count per gene" # otherwise problems with zero.replacement in next step
colnames(interactions.per.Mm.gene)[2] <- "Count per gene"

## replace zero with NA for easier exclusion from histogram
interactions.per.Ef.gene[interactions.per.Ef.gene == 0] <- NA
interactions.per.Mm.gene[interactions.per.Mm.gene == 0] <- NA


## Data prep for plotting
## efalci
df.ef <- na.omit(interactions.per.Ef.gene)
df.ef <- df.ef[with(df.ef, order(df.ef$`Count per gene`)), ]
df.ef$rownums <- nrow(df.ef):1

## mouse
df.mm <- na.omit(interactions.per.Mm.gene)
df.mm <- df.mm[with(df.mm, order(df.mm$`Count per gene`)), ]
df.mm$rownums <- nrow(df.mm):1

## plot mouse and parasite densities together
ef.plot <- ggplot(df.ef, aes(x=df.ef$rownums, y=df.ef$`Count per gene`)) +
  geom_bar(stat = "identity", size = 1, colour= "#e19e2e") +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Count per gene") +
  scale_x_discrete("Genes") +
  ggtitle("Interactions per parasite gene")
ggsave(file = "Supplement/SI_Ef.interaction.distribution.svg", 
       height = 10, width = 12, plot = ef.plot)

mm.plot <- ggplot(df.mm, aes(x=df.mm$rownums, y=df.mm$`Count per gene`)) +
  geom_bar(stat="identity", size = 1, colour="dark green") +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Count per gene") +
  scale_x_discrete("Genes") +
  ggtitle("Interactions per mouse gene")
ggsave(file = "Supplement/SI_Mm.interaction.distribution.svg", 
       height = 10, width = 12, plot = mm.plot)
#dev.off()

### Take out all interacting genes AND how many interactions they have
## Eimeria and mouse
interacting.genes.Ef <- subset(interactions.per.Ef.gene, !is.na(interactions.per.Ef.gene$`Count per gene`))
interacting.genes.Mm <- subset(interactions.per.Mm.gene, !is.na(interactions.per.Mm.gene$`Count per gene`))

df.ef[head(sort(interacting.genes.Ef[,2], decreasing = TRUE), n = 10), ]

## or access specific numer of interactions by: df.ef[max(df.ef[,3]),]
## in this case the highest number

########################## Network visualization attempt ####################################
library(graph)
library(Rgraphviz)
edges <- list(a=list(edges=2:3),
              b=list(edges=2:3),
              c=list(edges=c(2,4)),
              d=list(edges=1))
g <- new("graphNEL", nodes=letters[1:4], edgeL=edges,
         edgemode="directed")
plot(g)
