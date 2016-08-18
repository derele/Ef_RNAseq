## qPCR data from Simone and Annica
library(scales)

setwd("~/Ef_RNAseq/")
qpcr.df <- as.data.frame(read.csv("data/NMRI1st2ndqPCR_frSS_EfctRef_long.csv"))

#q.ymin <- subset(qpcr.df, qpcr.df$gene %in% "Ef18S")[,7]-
#  subset(qpcr.df, qpcr.df$gene %in% "Ef18S")[,8] 

#q.ymax <- subset(qpcr.df, qpcr.df$gene %in% "Ef18S")[,7]+
#  subset(qpcr.df, qpcr.df$gene %in% "Ef18S")[,8] 

long.qpcr <- melt(qpcr.df, measure.vars = c("rep1_ct", "rep2_ct", "rep3_ct"))
## Normalise to highest ct among all replicates
long.qpcr$norm.value <- long.qpcr$value/max(long.qpcr$value)


## Calculate sd and means on normalised values
stats.qpcr <- summaryBy(norm.value ~ dpi + inf,
          data = long.qpcr, 
          FUN = list(mean, min, max, sd))
stats.qpcr <- merge(stats.qpcr, long.qpcr)

################ STATUS ####################
## No transformation of data needed
## plot on log2 scale: transformation on ggplot possible, but data is excluded 
## (due to neg values?): What do about neg values then? Are they "wrong"?

## Discussed Totta & Emanuel: EH will write now ~ 2 weeks
## he will focus more on upcoming grant application towards end of August
## Cool that we have so much qPCR data: find those genes in RNAseq data and 
## plot together

bars.qpcr <- ggplot(subset(stats.qpcr, stats.qpcr$gene %in% "Ef18S"), 
       aes(x = dpi, y = norm.value.mean, col = factor(inf))) +
  #geom_boxplot(aes(group = round_any(dpi, 2, floor))) +
  geom_errorbar(aes(ymin = norm.value.min,
                    ymax = norm.value.max,
                    width = 0.2)) +
  geom_point(aes(size = 2)) +
  scale_y_continuous("Ef18S log2 fold change (delta-delta-ct)") #, #trans=log2_trans()) + #, limits = c(1, 200000)) +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        #axis.line = element_line(colour = "black"),
        #axis.title.y = element_text(),
        axis.title.x = element_text("Day post infection"),
        legend.key = element_blank(),
        legend.text = element_blank())
  

bars.qpcr

