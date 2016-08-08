library(rtracklayer)

messed.gff <- import.gff("/SAN/Eimeria_Totta/reference_genomes/mm10_eimeria_merge/indexes_bowtie2/mm10_GRCm38_eimeriaHaberkorn.gtf")

messed.gff$gene_id <- gsub("_\\d+$", "", messed.gff$gene_id)
messed.gff$gene_name <- gsub("_\\d+$", "", messed.gff$gene_name)

good.gff <- messed.gff

export.gff(good.gff, "/SAN/Eimeria_Totta/reference_genomes/mm10_eimeria_merge/indexes_bowtie2/mm10_GRCm38_eimeriaHaberkorn_fixed.gtf")
