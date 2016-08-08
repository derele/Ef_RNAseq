library(Rsubread)
library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)

dirs <- list.dirs("/SAN/Eimeria_Totta/tophat_March_c", recursive=FALSE)

bams <- sapply(dirs, list.files, "accepted_hits.bam", full.names=TRUE)

FC <- featureCounts(bams,
                    annot.ext="/SAN/Eimeria_Totta/reference_genomes/mm10_eimeria_merge/indexes_bowtie2/mm10_GRCm38_eimeriaHaberkorn_fixed.gtf",
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE,
                    nthreads=20)

## the experimental meta-data as in the pData slot of all ballgown object
create.pdata <- function (featurecountsTarget){
    samples <- gsub(".*?tophat_March_c\\.(.*?_rep\\d).*", "\\1", 
                   featurecountsTarget)
    grouped <- sub("_rep\\d+$" , "", samples)
    infection <- ifelse(grepl("1stInf", grouped), "1st",
                        ifelse(grepl("2ndInf", grouped), "2nd", "NULLth"))
    rep <- sub(".*_(rep\\d+)$" , "\\1", samples) # (name here) is taken "\\1" here
    mouse.strain <- sub("^(.*?)_.*", "\\1", grouped)
    timepoint <- as.numeric(as.character(sub("^.*_(.*)dpi", "\\1", grouped)))
    ## variables to summarize other variables sith mutliple levels into two levels ##
    immune.status <- ifelse(mouse.strain %in% "Rag", "Rag", "competent")
    late.early <- ifelse(timepoint < 7 , "early", "late")
    batch <- c(3, 3, 3, 3, 3, 1, 1, 2, 2, 2, 0, 1, 1, 2,
               3, 3, 3, 2, 3, 3, 3, 1, 0, 1, 0, 3, 3, 3, 3, 3)
    seq.method <- c("hiseq", "hiseq", "hiseq", "hiseq", "hiseq",
               "GAII", "GAII", "GAII", "GAII", "GAII", "GAII", "GAII", "GAII", "GAII",
               "hiseq", "hiseq", "hiseq",
               "GAII",
               "hiseq", "hiseq", "hiseq",
               "GAII", "GAII", "GAII", "GAII",
               "hiseq", "hiseq", "hiseq", "hiseq", "hiseq")
    sample.seq <- paste(samples,seq.method)
    data.frame(samples, grouped, infection, rep, mouse.strain,
               timepoint, immune.status, late.early, batch, seq.method, sample.seq)
    
}

pData <- create.pdata(FC$targets)

FC$targets <- gsub(".*?tophat_March_c\\.(.*?_rep\\d).*", "\\1", FC$targets)
colnames(FC$counts) <- gsub(".*?tophat_March_c\\.(.*?_rep\\d).*", "\\1", colnames(FC$counts))

write.table(FC$counts, "output_data/RC_All_genes.csv", sep=",")

write.table(pData, "output_data/Sample_pData.csv", sep=",")

