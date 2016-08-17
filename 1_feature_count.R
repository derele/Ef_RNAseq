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


colnames(FC$counts) <- gsub(".*?tophat_March_c\\.(.*?_rep\\d).*", "\\1", colnames(FC$counts))

## The phenotype information stored in fastq and bam file/folder names
## is inconsistent in that the following samples lack challenge/first
## infection information compared to e.g. NMRI_2ndInf_0dpi_rep1

## They are 1st_Inf and 2nd_Inf controls a lacking for these
## conditions.

colnames(FC$counts)[colnames(FC$counts)%in%"Rag_0dpi_rep1"] <-
    "Rag_1stInf_0dpi_rep1"
colnames(FC$counts)[colnames(FC$counts)%in%"Rag_0dpi_rep2"] <-
    "Rag_1stInf_0dpi_rep2"
colnames(FC$counts)[colnames(FC$counts)%in%"C57BL6_0dpi_rep1"] <-
    "C57BL6_1stInf_0dpi_rep1"
colnames(FC$counts)[colnames(FC$counts)%in%"C57BL6_0dpi_rep2"] <-
    "C57BL6_1stInf_0dpi_rep2"

write.table(FC$counts, "output_data/RC_All_genes.csv", sep=",")


