## Summarize and create plots or tables for data
## E.g. raw counts, percetage Eimeria reads... etc

library(scales)
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
#colnames(tableEf)[,1] <- "Read.number" - doesnt work?
#tableEf$Samples <- colnames(Ef.RC[[3]])
#print.xtable(tableEf, type = "latex", file = "tableEf-reads.tex")
tableEf$Samples <- rownames(tableEf)
## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]
row.names(tableEfsmall) <- seq(1:18)

## ADD EIMERIA % to EIMERIAsmall table
tableEfsmall$percent.Ef.reads <- (tableEfsmall[,1]/tableMsmall[,1])*100

png("figures/ef.read.numbers.png", width = 40, hight = 40, dpi = 300)
ef.read.number <- ggplot()+
	geom_point(data = tableEf, aes(x = Samples, y = as.numeric(colSums(Ef.RC[[3]])))) +
	scale_y_log10(labels = comma) +
	ylab("Number of Eimeria reads") +
	ggtitle("Sequences per samples for Eimeria") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ef.read.numbers.png", ef.read.number, path = "figures")#, width = 40, hight = 40, dpi = 300)
#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))







