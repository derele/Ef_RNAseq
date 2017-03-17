library(ggplot2)
library(reshape)
library(gridExtra)
library(xtable)
library(RSvgDevice)

## a) decide for samples needed to be excluede because of uncertain
## infection status

## b) create Table one summarizing metadata (previously pData) and
## basic sequencing statistics read counts

## Raw counst for tenes
All.RC <- as.matrix(read.table("output_data/RC_All_genes.csv", sep=","))

create.pdata <- function (featurecountsTarget){
    sample <- featurecountsTarget
    grouped <- sub("_rep\\d+$" , "", sample)
    challenged <- unlist(lapply(strsplit(sample, "_"), "[[", 2))
    rep <- sub(".*_(rep\\d+)$" , "\\1", sample) # (name here) is taken "\\1" here
    mouse.strain <- sub("^(.*?)_.*", "\\1", grouped)
    dpi <- gsub(".*?_(.*?)$", "\\1", grouped)
    ## below information that can't be extracted out ouf the truncated
    ## mapping file name
    batch <- c(3, 3, 3, 3, 3, 1, 1, 2, 2, 2, 0, 1, 1, 2,
               3, 3, 3, 2, 3, 3, 3, 1, 0, 1, 0, 3, 3, 3, 3, 3)
    seq.method <- c("hiseq", "hiseq", "hiseq", "hiseq", "hiseq",
                    "GAII", "GAII", "GAII", "GAII", "GAII", "GAII",
                    "GAII", "GAII", "GAII", "hiseq", "hiseq", "hiseq",
                    "GAII", "hiseq", "hiseq", "hiseq", "GAII", "GAII",
                    "GAII", "GAII","hiseq", "hiseq", "hiseq", "hiseq",
                    "hiseq")
    data.frame(sample, grouped, challenged, rep, mouse.strain,
               dpi, batch, seq.method)
    
}

RC.table <- create.pdata(colnames(All.RC))

## get the total read count.
## see README_fastq_linecount.txt on how this file was produced    
rawfastqRC <- readLines("/SAN/Eimeria_Totta/fastq_linecounts.txt")

fastqRC <- rawfastqRC[seq(1, length(rawfastqRC), by=2)+1]
fastqRC <- as.numeric(as.character(fastqRC))/4

names(fastqRC) <- rawfastqRC[seq(1, length(rawfastqRC), by=2)]
names(fastqRC) <- gsub("^.*?\\/(.*rep\\d).*", "\\1", names(fastqRC) )

## The nasty file renaming here as well
names(fastqRC)[names(fastqRC)%in%"Rag_0dpi_rep1"] <-
    "Rag_1stInf_0dpi_rep1"
names(fastqRC)[names(fastqRC)%in%"Rag_0dpi_rep2"] <-
    "Rag_1stInf_0dpi_rep2"
names(fastqRC)[names(fastqRC)%in%"C57BL6_0dpi_rep1"] <-
    "C57BL6_1stInf_0dpi_rep1"
names(fastqRC)[names(fastqRC)%in%"C57BL6_0dpi_rep2"] <-
    "C57BL6_1stInf_0dpi_rep2"

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

## plotting parasite percentage of reads

#pdf("figures/Figure1b_Ef.percentage.pdf")
#png("figures/Figure1b_Ef.percentage.png")

## If plot title is wanted:
#my.title <- expression(paste("Fraction of total ", italic("Eimeria"), " sequences per sample"))
my.ylab <- expression(paste("Percentage of reads mapping to", italic("Eimeria"), " genome"))

pdf("figures/Figure1d_Ef.NMRI.percentage.pdf", height = 14, width = 25)
ggplot(RC.table[RC.table$mouse.strain %in% "NMRI", ], 
       aes(x = dpi, y = p.Ef.reads,
                     color=challenged, shape=strain)) +
    geom_point(size=6)+
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides="lr") +
    theme_bw(32) +
    theme(axis.text.x = element_text(hjust = 0.5), #, size = 32),
          axis.title.x = element_text(size = 32),
          #axis.text.y = element_text(size = 32),
          axis.title.y = element_text(size = 32),
          title = element_text(size = 40)) +
    #theme_bw() +
    theme(legend.key = element_blank()) +
    ylab(label = my.ylab) +
    xlab(label = "Day post infection")
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

pdf("Supplement/FigureS2_distributionsMm.pdf", width = 27, height = 21)
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
               
pdf("Supplement/FigureS2_distributionsEf.pdf", width = 27, height = 21)
do.call(grid.arrange, c(density.plots, list(nrow=2)))
dev.off()


## a cut-off of e.g 100 is enough for Eimeria

## but we need to exclude samples "NMRI_2ndInf_3dpi_rep1" and
## "NMRI_2ndInf_5dpi_rep2" both samples are below the number of reads
## (and percentage) of the "NMRI_1stInf_0dpi_rep1" sample which is
## excluded for uncertain infection status. See:
RC.table <- RC.table[order(RC.table$c.Ef.reads),]

write.table(RC.table, "output_data/Sample_pData.csv", sep=",")


## I tink we don't need dpi, and challenged, this is in the Sample
## name and can be explained in a subtext
table.cleaned <- RC.table[, c("sample", "seq.method", "batch",
                              "c.seq.reads", 
                              "c.Mm.reads", "c.Ef.reads",
                              "p.Ef.reads", "c.Ef.genes")]

## Pretty names for tex
names(table.cleaned) <- sub("^seq.method$", "Sequencing method", names(table.cleaned))
names(table.cleaned) <- sub("^sample$", "Sample", names(table.cleaned))
names(table.cleaned) <- sub("^c.seq.reads$", "total reads", names(table.cleaned))
names(table.cleaned) <- sub("^c.Mm.reads$", "reads mapping Mouse", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.reads$", "reads mapping E. falciformis", names(table.cleaned))
names(table.cleaned) <- sub("^p.Ef.reads$", "Percentage E. falciformis", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.genes$", "# E. falciformis genes", names(table.cleaned))

## EXPORT to Latex format
table.tex <- xtable(table.cleaned, align = c("l", "l", "l", "l", "l", "l", "l", "l", "l"),
                    digits = c(0, 0, 0, 0, 0, 0, 0, 4, 0))
## for LaTeX output
#print(table.tex, type = "latex", file = "tables/Table1_ReadCounts.tex", include.rownames = F,
#      format.args = list(big.mark = ",", decimal.mark = "."))

## for HTML output and usage in e.g. Word, Open Office...
print(table.tex, type = "html", file = "tables/Table1_ReadCounts.html", include.rownames = F,
      format.args = list(big.mark = ",", decimal.mark = "."))



