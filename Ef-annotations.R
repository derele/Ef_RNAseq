# Annotations Eimeria falciformis

go_Ef <- read.csv("~/12864_2014_6777_GO-annot.csv")


############ below line does not worki#################
#flat.GO <- lapply(SY2GO, function (x) do.call("rbind", x))

#GO.df <- as.data.frame(do.call("rbind", flat.GO))
#GO.df$gene_id <- rep(names(flat.GO), times = unlist(lapply(flat.GO, nrow)))

gene2GO_Ef <- by(go_Ef, go_Ef$gene, function(x) as.character(x$go))

## subset for only those genes known that have been tested to create
## the suitable universe to search enrichment against
all.genes.E <- unique(select.from.stats.results(stat_results_time, Ef.bg,
                                              qval = 2))

##### NOTES ################
# What do I need from the 12....csv files in the home folder?
# Are the GO:ids enough for the enrichment tests
# Is the format correct now?
# Where do the _1 _2 come from and how do we deal with them?
# Why don't the .pdf rep-pair-comparisons work? (png mouse is fine) 
# How to implement the FCS(?) method

#### Stuff to consider - but also potential
#### time thieves:
# Plot No. of zero read genes VS total?
# Compare distributions in microarray VS RNAseq?
#### SO:
# Play with the 100-threshold differently
# for mouse and Eimeria.
# If necessary: invest more in finding comments
# about edgeR and bimodal distributions/the neg. binom distr. 

## adjust this for Ef
#all.genes <- get.annotation.for.xloc(all.genes)[[2]]
#gene2GO <- gene2GO[names(gene2GO)%in%all.genes]

