## qPCR data from Simone and Annica
library(scales)
library(reshape)
library(svglite)

setwd("~/Ef_RNAseq/")
## Ef18S qPCR data for different timepoints, 1st and 2nd infection and 3 replicates.
qpcr.df <- as.data.frame(read.csv("data/NMRI1st2ndqPCR_frSS_EfctRef_long.csv"))
long.qpcr <- melt(qpcr.df, measure.vars = c("rep1_ct", "rep2_ct", "rep3_ct"))
## Normalise to highest ct among all replicates
long.qpcr$norm.value <- (long.qpcr$value - max(long.qpcr$value))*-1


## Calculate sd and means on normalised values
stats.qpcr <- summaryBy(norm.value ~ dpi + inf,
          data = long.qpcr, 
          FUN = list(mean, min, max, sd))
stats.qpcr <- merge(stats.qpcr, long.qpcr)

################ STATUS ####################
## No transformation of data needed:
## plot on log2 scale: transformation on ggplot possible, but data is excluded 
## (due to neg values?): Negative values = more of target gene than of reference(?)

## Discussed Totta & Emanuel: EH will write now ~ 2 weeks
## he will focus more on upcoming grant application towards end of August
## Cool that we have so much qPCR data: find those genes in RNAseq data and 
## plot together

qpcr18S <- ggplot(subset(stats.qpcr, stats.qpcr$gene %in% "Ef18S"), 
       aes(x = dpi, y = norm.value.mean, col = factor(inf))) +
  geom_errorbar(aes(ymin = norm.value.min,
                    ymax = norm.value.max,
                    width = 0.2)) +
  geom_point(aes(size = 2)) +
  ##  (delta-delta-ct): put in figure legend
  scale_y_continuous("Eimeria 18S fold change", labels = math_format(2^.x), breaks = seq(0,2^30, 2)) +
  scale_x_continuous("Day post infection", breaks = seq(0,8,1)) +
  theme_bw(20) +
  theme(#legend.title = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggsave("figures/Figure1c_qPCR_Ef18S.pdf", plot = qpcr18S)
qpcr18S


