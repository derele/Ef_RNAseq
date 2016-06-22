#### WRONG
## Something is wrong here: use "A_data_curation.R" instead!!!!!!! 

#controlsc## Summarize and create plots or tables for data
## E.g. raw counts, percetage Eimeria reads... etc

library(scales)
library(xtable)
library(stringr)

#<<<<<<< HEAD
####################################################################
# TABLE WITH RAW DATA - OUTPUT IN LATEX
####################################################################
#=======
## Summary table of transcript counts per sample
## MOUSE table
tableM <- as.data.frame(colSums(Mm.RC[[3]]))
#colnames(tableM)[,2] <- "Read.number"
## MOUSE w/o day 0 
tableM$Samples <- rownames(tableM)

tableMsmall <- tableM[!(str_detect(rownames(tableM),"^.*0dpi.*$")), ]
row.names(tableMsmall) <- seq(1:17)
#>>>>>>> Totta_working

## table MOUSE and EIMERIA raw counts and percentage parasite reads 
table.rc <- as.data.frame(colSums(All.RC[[3]][str_detect(rownames(All.RC[[3]]), "^EfaB_.*$"), ]))
#tableEf <- as.data.frame(colSums(Ef.RC[[3]]))
colnames(table.rc)[1] <- "EimeriaTranscripts"
table.rc$MusTranscripts <- colSums(All.RC[[3]][str_detect(rownames(All.RC[[3]]), "^ENSMUS.*$"), ])
table.rc$Samples <- rownames(table.rc)
row.names(table.rc) <- seq(1:30)
table.rc$PercentEimeria <- format((table.rc[,1]/colSums(All.RC[[3]]))*100, digits=3, scientific = F)

#<<<<<<< HEAD
# sort and make rownames pretty for LaTeX use
sorted.table.rc <- table.rc[order(table.rc$PercentEimeria, decreasing = T), ]
sorted.table.rc <- sorted.table.rc[, c(3,4,2,1)]

## EXPORT to Latex format
rawcount.table <- xtable(sorted.table.rc, align = c("l", "l", "l", "l", "l"), digits = 3)
print(rawcount.table, type = "latex", file = "figures/table_rc.tex", include.rownames = F)

## Same as above but with order of magnitutde in numbers
rawcount.table2 <- xtable(sorted.table.rc, align = c("l", "l", "l", "l", "l"), digits = -3)
print(rawcount.table2, type = "latex", file = "figures/table_rc2.tex", include.rownames = F)


## add timepoint here for plots - but not above in table
table.rc$Timepoint <- pData(Ef.bg)$timepoint 

## plotting parasite percentage of reads
ef.percentage <- ggplot() +
  geom_point(data = table.rc, aes(x = Samples, y = PercentEimeria)) + #, color = pData(Ef.bg)[, 7])) +
#=======
## ADD EIMERIA % to EIMERIAsmall table
as.data.frame(colSums(Ef.RC[[3]]))
#print.xtable(tableEf, type = "latex", file = "tableEf-reads.tex")
tableEf$Samples <- rownames(tableEf)

## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]
row.names(tableEfsmall) <- seq(1:17)
tableEfsmall$mouseReads <- cbind(tableMsmall[,1])

## add oocysts and sporozoites
oospor <- tableEf[16:17, ]
oospor$mouseReads <- c(0, 0)
#oospor$Percent.Eimeria.transcripts <- c(99.9, 99.9)
#colnames(oospor)[1] <- "Transcripts"
row.names(oospor) <- seq(19:20)
tableEfsmall <- rbind(tableEfsmall, oospor)
## rearrange columns
tableEfsmall <- tableEfsmall[, c(2,3,1)]

## ADD EIMERIA % to EIMERIAsmall table
tableEfsmall$PercentEimeria <-(tableEfsmall[,3]/(tableEfsmall[,2]+tableEfsmall[,3]))*100
#tableEfsmall$percent.Ef.reads <- 
#      (tableEfsmall[,1]/(tableMsmall[,1]+tableEfsmall[,1]))*100
colnames(tableEfsmall) <- c("Samples", "MusTranscripts", "EimeriaTranscripts", "PercentEimeria")
# sort and make rownames pretty for LaTeX use
sortedEfsmall <- tableEfsmall[order(tableEfsmall$PercentEimeria, decreasing=T), ]
rownames(sortedEfsmall) <- seq(1:20)

#tableEfsmall <- format(tableEfsmall, scientific =F, digits=1, justify = "centre")
efsmall <- xtable(sortedEfsmall, align = c("c", "c", "c", "c", "c"), digits = 0)
print(efsmall, type = "latex", file = "figures/tableEfsmall-2.tex")
## put percentage into table (plot not ideal due to outliers)

## add timepoint here for plots - but not above in table
tableEfsmall$Timepoint <- c(5, 5, 5, 3, 3, 5, 5, 7, 7, 3, 3, 5, 5, 7, 7, 5, 5, 5, NA, NA)

## plotting parasite percentage of reads
ef.percentage <- ggplot() +
  geom_point(data = tableEfsmall, aes(x = Samples, y = PercentEimeria)) + #, color = pData(Ef.bg)[, 7])) +
#>>>>>>> Totta_working
  ylim(c(0, 10))+
  ylab("Percentage of reads mapping to Eimeria genome") +
  ggtitle("Eimeria fraction of total sequences per sample - 2 samples rm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.percentage-cutoff10perc.png", ef.percentage, path = "figures", width = 40)
dev.off()
	
## plot number of parasite reads
ef.read.number <- ggplot()+
#<<<<<<< HEAD
  geom_point(data = table.rc, aes(x = Samples, y = as.numeric(EimeriaTranscripts), color = factor(Timepoint))) +
#=======
  geom_point(data = tableEfsmall, aes(x = Samples, y = as.numeric(EimeriaTranscripts), color = factor(Timepoint))) +
#>>>>>>> Totta_working
  scale_y_log10(labels = comma) +
  ylab("Number of Eimeria reads") +
  ggtitle("Sequences per samples for Eimeria - immune-status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.read.numbers-immune-status.png", ef.read.number, path = "figures", width = 40)#, hight = 40, dpi = 300)
dev.off()
#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))

