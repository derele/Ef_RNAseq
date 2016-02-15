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
ts.data.frame(colSums(Ef.RC[[3]]))
#print.xtable(tableEf, type = "latex", file = "tableEf-reads.tex")
tableEf$Samples <- rownames(tableEf)

## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]
row.names(tableEfsmall) <- seq(1:18)

## ADD EIMERIA % to EIMERIAsmall table
tableEfsmall$percent.Ef.reads <- (tableEfsmall[,1]/tableMsmall[,1])*100
colnames(tableEfsmall) <- c("Transcripts", "Samples", "Percent.Eimeria.transcripts")

#tableEfsmall <- format(tableEfsmall, scientific =F, digits=1, justify = "centre")
efsmall <- xtable(tableEfsmall, align = c("c", "c", "c", "c"), digits = 1)
print(efsmall, type = "latex", file = "figures/tableEfsmall.tex")
## put percentage into table (plot not ideal due to outliers)



## plotting parasite percentage of reads
ef.percentage <- ggplot() +
  geom_point(data = tableEfsmall, aes(x = Samples, y = Percent.Eimeria.transcripts)) + #, color = pData(Ef.bg)[, 7])) +
  ylim(c(0, 10))+
  ylab("Percentage of reads mapping to Eimeria genome") +
  ggtitle("Eimeria fraction of total sequences per sample - 2 samples rm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.percentage-cutoff10perc.png", ef.percentage, path = "figures", width = 40)
dev.off()
	
## plot number of parasite reads
ef.read.number <- ggplot()+
  geom_point(data = tableEf, aes(x = Samples, y = as.numeric(colSums(Ef.RC[[3]])), color = pData(Ef.bg)[ , 7] )) +
  scale_y_log10(labels = comma) +
  ylab("Number of Eimeria reads") +
  ggtitle("Sequences per samples for Eimeria - immune-status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.read.numbers-immune-status.png", ef.read.number, path = "figures", width = 40)#, hight = 40, dpi = 300)
dev.off()
#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))

