## all simple functions go here

TOGO.all.onto <- function (ontology, allgenes, gene.set, annot) {
  g <- factor(as.integer( allgenes%in%gene.set ))
  names(g) <- allgenes
  toGO <-  new("topGOdata", ontology = ontology, allGenes = g, annot = annFUN.gene2GO,
               gene2GO = annot)
  resultFis <- runTest(toGO, algorithm = "classic", statistic = "fisher")
  list(toGO, resultFis) ## returns a list first data then result
}


gene.table.topGO <- function(TOGO.list, pval=0.01){
    all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
    ##    all$fdr <- p.adjust(all$result1, method="BH")
    names(all)[names(all)%in%"result1"] <- "p.value"
    all$p.value <- gsub(" ?< ?", "", all$p.value)
    all$p.value <- as.numeric(all$p.value)
    all$adj.p <- round(p.adjust(all$p.value, method="BH"), 4)
    return(all[all$p.value<pval,])
}


get.genes.4.clusters <- function (clusters) {
  gsub("Efa\\|", "",   grep("Efa\\|",
                            unlist(clusters),
                            value=TRUE))
}

## adjusted functions from ballgown
#if(!require(GenomicRanges)) biocLite("GenomicRanges")

library(ballgown)
library(GenomicRanges) # needed for elementMetadata()

############################################################################################
#### change SUBSET function to allow logical input for cond argument
setMethod("subset", "ballgown", function(x, cond, genomesubset=TRUE){
  stopifnot(class(cond) == 'logical') # CHANGED TO LOGICAL
  # if you are subsetting by something in the genome (say, a chromosome):
  if(genomesubset){
    trans = subset(expr(x)$trans, cond)  # removed eval().... part
    thetx = trans$t_id
    
    inttmp = split(indexes(x)$i2t$i_id, indexes(x)$i2t$t_id)
    theint = as.numeric(unique(unlist(inttmp[names(inttmp) %in% thetx])))
    intron = subset(expr(x)$intron, i_id %in% theint)
    
    extmp = split(indexes(x)$e2t$e_id, indexes(x)$e2t$t_id)
    theex = as.numeric(unique(unlist(extmp[names(extmp) %in% thetx])))
    exon = subset(expr(x)$exon, e_id %in% theex)
    
    e2t = subset(indexes(x)$e2t, t_id %in% thetx)
    i2t = subset(indexes(x)$i2t, t_id %in% thetx)
    t2g = subset(indexes(x)$t2g, t_id %in% thetx)
    
    introngr = structure(x)$intron[elementMetadata(structure(x)$intron)$id 
                                   %in% theint]
    exongr = structure(x)$exon[elementMetadata(structure(x)$exon)$id 
                               %in% theex]
    grltxids = as.numeric(names(structure(x)$trans))
    transgrl = structure(x)$trans[grltxids %in% thetx]
    
    return(new("ballgown", expr=list(intron=intron, exon=exon, trans=trans),
               indexes=list(e2t=e2t, i2t=i2t, t2g=t2g, 
                            bamfiles=indexes(x)$bamfiles, pData=indexes(x)$pData), 
               structure=list(intron=introngr, exon=exongr, trans=transgrl), 
               dirs=dirs(x), mergedDate=mergedDate(x), meas=x@meas, RSEM=x@RSEM))
  }
  else{
    # If genomesubset is FALSE:
    # you're doing a phenotype subset
    # structure, some indexes, dirs, and mergedDate stay the same
    # change: data, indexes(pData), and indexes(bamfiles)
    ## pData
    newpd = subset(pData(x), cond) # removed eval().... part
    newpd = droplevels(newpd)
    ## bamfiles
    newsampnames = newpd[,1]
    rowIndsToKeep = which(pData(x)[,1] %in% newsampnames)
    if(!is.null(indexes(x)$bamfiles)){
      newbamfiles = indexes(x)$bamfiles[rowIndsToKeep]
    }else{
      newbamfiles = NULL
    }
    ## transcript data
    ## we changed the column selection proceedure completely ;-)...
    txKeepCols = c(rep(TRUE, times=10), rep(cond, each = 2))
    newtdat = texpr(x, 'all')[,txKeepCols] 
    ## exon data
    exKeepCols = c(rep(TRUE, times=5), rep(cond, each = 7))
    newedat = eexpr(x, 'all')[,exKeepCols] 
    ## intron data
    iKeepCols = c(rep(TRUE, times=5), rep(cond, each = 3))
    newidat = iexpr(x, 'all')[,iKeepCols] 
    return(new("ballgown", 
               expr=list(intron=newidat, exon=newedat, trans=newtdat), 
               indexes=list(e2t=indexes(x)$e2t, i2t=indexes(x)$i2t,
                            t2g=indexes(x)$t2g, bamfiles=newbamfiles, pData=newpd), 
               structure=list(intron=structure(x)$intron, exon=structure(x)$exon, 
                              trans=structure(x)$trans),
               dirs=dirs(x)[rowIndsToKeep], mergedDate=mergedDate(x), 
               meas=x@meas, RSEM=x@RSEM))
  }
} )

#################################################################################################
