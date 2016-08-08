## This is a script is 2nd the pipeline and needed for evaluation
## purpose

## Goals:

## EX- on IN-clusion of samples and genes ###################

## a) decide for samples needed to be excluede because of uncertain
## infection status

## b) decice for a filter on minimal expression values over all
## samples to exclude genes

## c) decide on samples needed to be excluded because of bad
## "unnormalizable" distribution of expression values


### Normalization

## d) decide on a general normalization strategy scaling the smaple
## throughput

## e) decide on a "removeal of unwanted variation" strategy

## f) decide on a potential addional "batch effect" removal strategy


## a) decide for samples needed to be excluede because of uncertain
## infection status

## The raw counts  required from upstream scripts:
## GENES
## All.RC <- as.matrix(read.table("output_data/RC_All_genes.csv", sep=","))

## TRANSCRIPTS
All.RC <- as.matrix(read.table("output_data/RC_All_genes.csv", sep=","))

## The phenotype information
All.pData <- read.table("output_data/Sample_pData.csv", sep=",")

library(ggplot2)
library(reshape)
library(gridExtra)
library(xtable)

## Overview for differential representation of host/parasite
RC.table <-  as.data.frame(cbind(sample=as.character(All.pData$samples),
                                 seq.method=as.character(All.pData$seq.method),
                                 batch=All.pData$batch))

## get the total read count.
## see README_fastq_linecount.txt on how this file was produced    
rawfastqRC <- readLines("/SAN/Eimeria_Totta/fastq_linecounts.txt")

fastqRC <- rawfastqRC[seq(1, length(rawfastqRC), by=2)+1]
fastqRC <- as.numeric(as.character(fastqRC))/4

names(fastqRC) <- rawfastqRC[seq(1, length(rawfastqRC), by=2)]
names(fastqRC) <- gsub("^.*?\\/(.*rep\\d).*", "\\1", names(fastqRC) )

RC.table <- merge(RC.table, fastqRC, by.x="sample", by.y=0)
names(RC.table)[ncol(RC.table)] <- "c.seq.reads"

## this is total reads mapped not total reads!!!
RC.table$c.mapping.counts <- colSums(All.RC)

## wtf more reads mapped than sequenced
RC.table$strain <- unlist(lapply(strsplit(as.character(RC.table$sample), "_"), "[[", 1))

RC.table$c.Mm.reads <- colSums(All.RC[grepl('^ENSMUS.*', rownames(All.RC)),])

RC.table$c.Ef.reads <- colSums(All.RC
                                [grepl('^EfaB.*', rownames(All.RC)),])

RC.table$p.Ef.reads <- round(RC.table$c.Ef.reads/RC.table$c.mapping.counts*100, 4)

RC.table$dpi <- gsub(".*?(\\ddpi).*", "\\1", RC.table$sample)

RC.table$challenged <- ifelse(grepl("2ndInf", RC.table$sample), "Callenge",
                       ifelse(grepl("0dpi", RC.table$sample), "None",
                       ifelse(grepl("1stInf", RC.table$sample), "First",
                              "environmental")))

RC.table$dpi[RC.table$challenged%in%"environmental"] <- "environmental"

## plotting parasite percentage of reads

pdf("figures/Figure1b_Ef.percentage.pdf")
ggplot(RC.table, aes(x = dpi, y = p.Ef.reads,
                     color=challenged, shape=strain)) +
    geom_jitter(size=4, width=0.35, alpha=0.6)+
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    annotation_logticks(sides="lr") +
    ylab(label="Percentage of reads mapping to Eimeria genome") +
    xlab(label="Time after infection") +
    ggtitle("Eimeria fraction of total sequences per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()
dev.off()

## taking a minimum of 5 reads as evidence of expression
RC.table$c.Ef.genes <- colSums(All.RC
                               [grepl('^EfaB.*', rownames(All.RC)),]>5)

## -> NMRI_1stInf_0dpi_rep1 has to be excluded from the mouse and
## parasite dataset because of uncretain infection status


## b) decice for a filter on minimal expression values over all samples

## For mouse
keep.val <- c(0, 100, 1000, 2000, 3000, 5000)
names(keep.val) <- keep.val

r.c.s.l <-
    lapply(keep.val,
           function (x){
               ## add tiny number to be able to plot transcripts with zero reads mapping
               Mm.RCtemp <- All.RC[grepl('^ENSMUS.*', rownames(All.RC)),]+ 0.1 
               kept.df <- Mm.RCtemp[rowSums(Mm.RCtemp)>x,]       
               kept.df <- melt(kept.df)
               ## remove ooc and sporo samples
               kept.df <- kept.df[grep("^.*_(oocysts|sporozoites)_.*$",
                                       kept.df$X2, invert = T),] 
               return(kept.df)
           })

#####################################################
density.plots <-
    lapply(seq_along(r.c.s.l),
           function(i){
               ggplot(r.c.s.l[[i]],
                      aes(value, ..density..)) + 
                   stat_density(geom="line") +
                   facet_wrap(~X2)+ 	# makes gray box on top of each plot
		   scale_x_log10("Read counts (log10)",
			labels = scales::trans_format("log10",
                                                      scales::math_format(10^.x))) +
	   	   theme(axis.text = element_text(size = 16),
			 axis.line = element_line(colour = "black", size=2),
			 plot.title = element_text(size = 20),
			 axis.title.x = element_text(size = 20),
			 axis.title.y = element_text(size = 20),
			 panel.background = element_blank(),
			 panel.grid.major = element_blank(),
			 panel.grid.minor = element_blank()) +
                   ggtitle(paste("Cutoff =",
                               names(r.c.s.l)[[i]]))
           })

pdf("Supplement/distributionsMm.pdf", width = 27, height = 21)
do.call(grid.arrange, c(density.plots, list(nrow=2)))
dev.off()

## For mouse we decide on a minimal value of 3000 reads across all
## samples to include a gene

## For Eimeria we need much less
keep.val <- keep.val/10
names(keep.val) <- keep.val

r.c.s.l <-
    lapply(keep.val,
           function (x){
               Ef.RCtemp <- All.RC[grepl('^EfaB.*', rownames(All.RC)),] + 0.1 # add tiny number to be able to plot transcripts with zero reads mapping	
               kept.df <- Ef.RCtemp[rowSums(Ef.RCtemp)>x,]       
               kept.df <- melt(kept.df)
               kept.df <- kept.df[grep("^.*_0dpi_.*$", kept.df$X2, invert = T),] # remove day 0 samples
               return(kept.df)
           })

density.plots <-
    lapply(seq_along(r.c.s.l),
           function(i){
               ggplot(r.c.s.l[[i]],
                      aes(value, ..density..)) + 
                   stat_density(geom="line") +
                   facet_wrap(~X2)+ 	# makes gray box on top of each plot
                   scale_x_log10("Read counts (log10)",
                                 labels = scales::trans_format("log10",
                                                               scales::math_format(10^.x))) +
	   	   theme(axis.text = element_text(size = 16),
			 axis.line = element_line(colour = "black", size=2),
			 plot.title = element_text(size = 20),
			 axis.title.x = element_text(size = 20),
			 axis.title.y = element_text(size = 20),
			 panel.background = element_blank(),
			 panel.grid.major = element_blank(),
			 panel.grid.minor = element_blank()) +
                   ggtitle(paste("Cutoff =",
                               names(r.c.s.l)[[i]]))
           })
               
pdf("Supplement/distributionsEf.pdf", width = 27, height = 21)
do.call(grid.arrange, c(density.plots, list(nrow=2)))
dev.off()


## a cut-off of e.g 100 is enough for Eimeria

## but we need to exclude samples "NMRI_2ndInf_3dpi_rep1" and
## "NMRI_2ndInf_5dpi_rep2" both samples are below the number of reads
## (and percentage) of the "NMRI_1stInf_0dpi_rep1" sample which is
## excluded for uncertain infection status. See:
RC.table <- RC.table[order(RC.table$c.Ef.reads),]

table.cleaned <- RC.table[, c("sample", "seq.method", "batch",
                              "c.Mm.reads", "c.Ef.reads",
                              "p.Ef.reads", "dpi", "challenged",
                              "c.Ef.genes")]

## Pretty names for tex
names(table.cleaned) <- sub("^seq.method$", "Sequencing method", names(table.cleaned))
names(table.cleaned) <- sub("^sample$", "Sample", names(table.cleaned))
names(table.cleaned) <- sub("^c.Mm.reads$", "read mappings Mouse", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.reads$", "read mappings E. falciformis", names(table.cleaned))
names(table.cleaned) <- sub("^p.Ef.reads$", "Percentage E. falciformis", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.genes$", "# E. falciformis genes", names(table.cleaned))

## EXPORT to Latex format
table.tex <- xtable(table.cleaned, align = c("l", "l", "l", "l", "l", "l", "l", "l", "l", "l"), digits = 3)
print(table.tex, type = "latex", file = "tables/Table1_ReadCounts.tex", include.rownames = F)



