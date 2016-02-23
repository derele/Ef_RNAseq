### ballgown for data import from the cufflinks/tablemaker pipeline

library(ballgown)
setwd("~/Ef_RNAseq")

## IMPORT DATA USING BALLGOWN IMPORT
if(!exists("All.bg")) {
    All.bg = ballgown(dataDir="/data/Eimeria_Totta/tablemaker_060815/",
        samplePattern = "*")
}

## some fundtions for RNAseq data, mainly a new subsetting function
## for ballgown-objects
source("functionsEfRNAseq.R")

########################

## the experimental meta-data to be stored in the pData slot of the object
add.pdata <- function (ballgown.obj){
    samples <- sub("^tablemakertophat_(.*)(_hiseq|_forw.fastq_paired)$",
                   "\\1", sampleNames(ballgown.obj))
    grouped <- sub("_rep\\d+$" , "", samples)
    infection <- ifelse(grepl("1stInf", grouped), "1st",
                        ifelse(grepl("2ndInf", grouped), "2nd", "NULLth"))
    rep <- sub(".*_(rep\\d+)$" , "\\1", samples) # (name here) is taken "\\1" here
    mouse.strain <- sub("^(.*?)_.*", "\\1", grouped)
    timepoint <- as.numeric(as.character(sub("^.*_(.*)dpi", "\\1", grouped)))
    ## variables to summarize other variables sith mutliple levels into two levels ##
    immune.status <- ifelse(mouse.strain %in% "Rag", "Rag", "competent")
    late.early <- ifelse(timepoint < 7 , "early", "late")
    batch <- c(3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 1,
               1, 2, 3, 3, 3, 2, 3, 3, 3, 1, 1, 3, 3, 3, 3, 3)
    seq.method <- c("hiseq", "hiseq", "hiseq", "hiseq", "hiseq",
                    "GAII", "GAII", "GAII", "GAII", "GAII", "GAII", "GAII", "GAII",
                    "hiseq", "hiseq", "hiseq", "GAII", "hiseq", "hiseq", "hiseq",
                    "GAII", "GAII", "hiseq", "hiseq", "hiseq", "hiseq", "hiseq")
    sample.seq <- paste(samples,seq.method)
    data.frame(samples, grouped, infection, rep, mouse.strain,
               timepoint, immune.status, late.early, batch, seq.method, sample.seq)
    
}

pData(All.bg) <- add.pdata(All.bg)

########################################
## Create mouse only object and remove sporozoite and oocyst samples
## from mouse data
Mm.bg <- subset(All.bg, grepl('^XLOC.*', geneIDs(rnaseq.bg)),
                genomesubset=TRUE)
Mm.bg <- subset(Mm.bg, !is.na(pData(Mm.bg)$timepoint),
                genomesubset=FALSE) # 

##################################################
## Create Eimeria only object and remove day zero
##################################################
Ef.bg <- subset(All.bg,
                grepl('^EfaB.*', geneIDs(All.bg)),
                genomesubset=TRUE)
Ef.bg <- subset(Ef.bg,
                pData(Ef.bg)$timepoint != 0 | is.na(pData(Ef.bg)$timepoint),
                genomesubset=FALSE)

## A function to et raw coungs for exons, transcripts and genes out of
## the ballgown objects

raw.counts.4.bg <- function(bg){
    exon.raw.count <- eexpr(bg, "ucount")
    ## the linkage data
    e2t <- bg@indexes$e2t
    t2g <- bg@indexes$t2g 
    e2t2g <- merge(e2t, t2g)
    all.count <- merge(e2t2g, exon.raw.count,
                       by.x = "e_id", by.y = 0)
    ## sum up for transcripts
    transcript.raw.count <-
        do.call("rbind",
                by(all.count, all.count$t_id,
                   function(x) colSums(x[, 4:ncol(x)])))
    ## sum up for genes
    gene.raw.count <-
        do.call("rbind",
                by(all.count, all.count$g_id,
                   function(x) colSums(x[, 4:ncol(x)])))
    colnames(transcript.raw.count) <-
        as.character(pData(bg)$samples)
    colnames(gene.raw.count) <-
        as.character(pData(bg)$samples)
    return(list(exon.raw.count,
                transcript.raw.count,
                gene.raw.count))
}

All.RC <- raw.counts.4.bg(All.bg)
Mm.RC <- raw.counts.4.bg(Mm.bg)
Ef.RC <- raw.counts.4.bg(Ef.bg)

