controlsc## Summarize and create plots or tables for data
## E.g. raw counts, percetage Eimeria reads... etc

library(scales)
library(xtable)
library(stringr)

####################################################################
# TABLE WITH RAW DATA - OUTPUT IN LATEX
####################################################################

## table MOUSE and EIMERIA raw counts and percentage parasite reads 
table.rc <- as.data.frame(colSums(All.RC[[3]][str_detect(rownames(All.RC[[3]]), "^EfaB_.*$"), ]))
#tableEf <- as.data.frame(colSums(Ef.RC[[3]]))
colnames(table.rc)[1] <- "EimeriaTranscripts"
table.rc$MusTranscripts <- colSums(All.RC[[3]][str_detect(rownames(All.RC[[3]]), "^XLOC_.*$"), ])
table.rc$Samples <- rownames(table.rc)
row.names(table.rc) <- seq(1:27)
table.rc$PercentEimeria <- format((table.rc[,1]/colSums(All.RC[[3]]))*100, digits=3, scientific = F)

# sort and make rownames pretty for LaTeX use
sorted.table.rc <- table.rc[order(table.rc$PercentEimeria, decreasing = T), ]
sorted.table.rc <- sorted.table.rc[, c(3,4,2,1)]

## EXPORT to Latex format
rawcount.table <- xtable(sorted.table.rc, align = c("c", "c", "c", "c", "c"), digits = 3)
print(rawcount.table, type = "latex", file = "figures/table_rc.tex", include.rownames = F)

## add timepoint here for plots - but not above in table
table.rc$Timepoint <- c(5, 5, 5, 3, 3, 5, 5, 7, 7, 3, 3, 5, 5, 7, 7, 5, 5, 5, NA, NA)

## plotting parasite percentage of reads
ef.percentage <- ggplot() +
  geom_point(data = table.rc, aes(x = Samples, y = PercentEimeria)) + #, color = pData(Ef.bg)[, 7])) +
  ylim(c(0, 10))+
  ylab("Percentage of reads mapping to Eimeria genome") +
  ggtitle("Eimeria fraction of total sequences per sample - 2 samples rm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.percentage-cutoff10perc.png", ef.percentage, path = "figures", width = 40)
dev.off()
	
## plot number of parasite reads
ef.read.number <- ggplot()+
  geom_point(data = table.rc, aes(x = Samples, y = as.numeric(EimeriaTranscripts), color = factor(Timepoint))) +
  scale_y_log10(labels = comma) +
  ylab("Number of Eimeria reads") +
  ggtitle("Sequences per samples for Eimeria - immune-status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.read.numbers-immune-status.png", ef.read.number, path = "figures", width = 40)#, hight = 40, dpi = 300)
dev.off()
#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))

