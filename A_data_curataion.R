## This is a script which is not part of the pipeline but needed for
## evaluation purpose -> only "insight" gained from this is transfered
## to the numberd scripts which are part of the processing pipeline

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

##object required from upstream scripts:

## The raw counts:
if(!exists("All.RC")|!exists("Mm.RC")|!exists("Ef.RC")){
    source("1_ballgown_import.R")
}

library(ggplot2)
library(reshape)
library(gridExtra)

## Overview for differential representation of host/parasite
RC.table <-  as.data.frame(cbind(sample=as.character(pData(All.bg)$samples),
                                 seq.method=as.character(pData(All.bg)$seq.method),
                                 batch=pData(All.bg)$batch))

RC.table$c.total.reads <- colSums(All.RC[[3]])

RC.table$c.Mm.reads <- colSums(All.RC[[3]][grepl('^ENSMUS.*', rownames(All.RC[[3]])),])

RC.table$c.Ef.reads <- colSums(All.RC[[3]]
                                [grepl('^EfaB.*', rownames(All.RC[[3]])),])

RC.table$p.Ef.reads <- round(RC.table$c.Ef.reads/RC.table$c.total.reads*100, 4)

## plotting parasite percentage of reads
ef.percentage <- ggplot(RC.table, aes(x = sample, y = p.Ef.reads,
                                      color=seq.method)) +
    geom_point()+
        scale_y_log10()+
            ylab(label="Percentage of reads mapping to Eimeria genome") +
                xlab(label="Sample") +
                ggtitle("Eimeria fraction of total sequences per sample") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("ef.percentage.png", ef.percentage, path = "figures", width = 10)

## taking a minimum of 5 reads as evidence of expression
RC.table$c.Ef.genes <- colSums(All.RC[[3]]
                               [grepl('^EfaB.*', rownames(All.RC[[3]])),]>5)

## -> NMRI_1stInf_0dpi_rep1 has to be excluded from the mouse and
## parasite dataset because of uncretain infection status


## b) decice for a filter on minimal expression values over all samples

## For mouse
keep.val <- c(0, 100, 1000, 2000, 3000, 5000)
names(keep.val) <- keep.val

## Original function
#r.c.s.l <-
#    lapply(keep.val,
#           function (x){
#               kept.df <- Mm.RC[[3]][rowSums(Mm.RC[[3]])>x, ]
#               kept.df <- melt(kept.df)
#               return(kept.df)
#           })

## Totta's fix to have data plotted with all (not already removed) samples:
r.c.s.l <-
    lapply(keep.val,
           function (x){
		
		Mm.RCtemp <- All.RC[[3]][grepl('^ENSMUS.*', rownames(All.RC[[3]])),]+ 0.1 # add tiny number to be able to plot transcripts with zero reads mapping
        	kept.df <- Mm.RCtemp[rowSums(Mm.RCtemp)>x,]       
		kept.df <- melt(kept.df)
		kept.df <- kept.df[grep("^.*_(oocysts|sporozoites)_.*$", kept.df$X2, invert = T),] # remove ooc and sporo samples
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
			labels = trans_format("log10",
			math_format(10^.x))) +
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
               
pdf("figures/distributionsMm.pdf", width = 27, height = 21)
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
		Ef.RCtemp <- All.RC[[3]][grepl('^EfaB.*', rownames(All.RC[[3]])),] + 0.1 # add tiny number to be able to plot transcripts with zero reads mapping	
        	kept.df <- Ef.RCtemp[rowSums(Ef.RCtemp)>x,]       
		kept.df <- melt(kept.df)
		kept.df <- kept.df[grep("^.*_0dpi_.*$", kept.df$X2, invert = T),] # remove day 0 samples
		return(kept.df)
           })
## OLD - remove after some days if everything works
#r.c.s.l <-
#    lapply(keep.val,
#           function (x){
#               kept.df <- Ef.RC[[3]][rowSums(Ef.RC[[3]])>x, ]
#               kept.df <- melt(kept.df)
#               return(kept.df)
#           })


density.plots <-
    lapply(seq_along(r.c.s.l),
           function(i){
               ggplot(r.c.s.l[[i]],
                   aes(value, ..density..)) + 
                   stat_density(geom="line") +
                   facet_wrap(~X2)+ 	# makes gray box on top of each plot
		   scale_x_log10("Read counts (log10)",
			labels = trans_format("log10",
			math_format(10^.x))) +
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
               
pdf("figures/distributionsEf.pdf", width = 27, height = 21)
do.call(grid.arrange, c(density.plots, list(nrow=2)))
dev.off()


foo <- melt(r.c.s.l)

pdf("figures/distributions_Ef_wrong.pdf")
ggplot(foo,
       aes(value, ..density.., color=X1)) + 
    stat_density(geom="line") +
        facet_wrap(~X2)+ 
            scale_x_log10("Read counts (log10)")



## a cut-off of e.g 100 is enough for Eimeria

## but we need to exclude samples "NMRI_2ndInf_3dpi_rep1" and
## "NMRI_2ndInf_5dpi_rep2" both samples are below the number of reads
## (and percentage) of the "NMRI_1stInf_0dpi_rep1" sample which is
## excluded for uncertain infection status. See:
RC.table[order(RC.table$c.Ef.reads),]

RC.table[order(RC.table$p.Ef.reads),]

table.cleaned <- RC.table[,c(1,5,6,7,8)]
## Pretty names for tex
names(table.cleaned) <- sub("^sample$", "Sample", names(table.cleaned))
names(table.cleaned) <- sub("^c.Mm.reads$", "Mouse transcripts", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.reads$", "E. falciformis transcripts", names(table.cleaned))
names(table.cleaned) <- sub("^p.Ef.reads$", "Percentage E. falciformis", names(table.cleaned))
names(table.cleaned) <- sub("^c.Ef.genes$", "#E. falciformis genes", names(table.cleaned))

table.cleaned <- table.cleaned[order(table.cleaned[,4], decreasing = T),]

## EXPORT to Latex format
table.tex <- xtable(table.cleaned, align = c("l", "l", "l", "l", "l", "l"), digits = 3)
print(table.tex, type = "latex", file = "figures/table_rc.tex", include.rownames = F)

table.cleaned2 <- table.cleaned[order(table.cleaned[,1]),]

## Same as above but only cols 1-3 and order of magnitutde in numbers
table.tex2 <- xtable(table.cleaned2[,1:3], align = c("l", "l", "l", "l"), digits = -3)
print(table.tex2, type = "latex", file = "figures/table_rc2.tex", include.rownames = F)


