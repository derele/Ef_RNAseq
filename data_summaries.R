## Summarize and create plots or tables for data
## E.g. raw counts, percetage Eimeria reads... etc

library(scales)
library(xtable)
library(stringr)

## Summary table of transcript counts per sample
## MOUSE table
tableM <- as.data.frame(colSums(mouse.RC[[3]]))
#colnames(tableM)[,2] <- "Read.number"
## MOUSE w/o day 0 
tableM$Samples <- rownames(tableM)

tableMsmall <- tableM[!(str_detect(rownames(tableM),"^.*0dpi.*$")), ]
row.names(tableMsmall) <- seq(1:18)

## EIMERIA table
tableEf <- as.data.frame(colSums(Ef.RC[[3]]))
tableEf$Samples <- rownames(tableEf)
## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]
row.names(tableEfsmall) <- seq(1:18)

## ADD EIMERIA % to EIMERIAsmall table
as.data.frame(colSums(Ef.RC[[3]]))
#print.xtable(tableEf, type = "latex", file = "tableEf-reads.tex")
tableEf$Samples <- rownames(tableEf)

## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]
row.names(tableEfsmall) <- seq(1:18)
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
  ylim(c(0, 10))+
  ylab("Percentage of reads mapping to Eimeria genome") +
  ggtitle("Eimeria fraction of total sequences per sample - 2 samples rm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.percentage-cutoff10perc.png", ef.percentage, path = "figures", width = 40)
dev.off()
	
## plot number of parasite reads
ef.read.number <- ggplot()+
  geom_point(data = tableEfsmall, aes(x = Samples, y = as.numeric(EimeriaTranscripts), color = factor(Timepoint))) +
  scale_y_log10(labels = comma) +
  ylab("Number of Eimeria reads") +
  ggtitle("Sequences per samples for Eimeria - immune-status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.read.numbers-immune-status.png", ef.read.number, path = "figures", width = 40)#, hight = 40, dpi = 300)
dev.off()
#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))

