## Reading and plotting Simones data on oocyst output and weigh change upon Eimeria infection
## of NMRI, C57BL/6 and Rag1-/- mice
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(reshape2)
library(plyr)
library(doBy)
library(stringr)
library(RSvgDevice)

setwd("~/Ef_RNAseq/")
phen.data <- read.csv("data/Oocysts_output_weight_SS_longdata.csv")

## Mouse IDs are not unique so it is impossible to track individuals FIX THIS!!
## here just a hack to take the first occurence of a ID at a time point
phen.data <- phen.data[!duplicated(phen.data[, c("Day_pi", "Mouse_ID", "Infection_No")]), ]

all.oocysts.line <-
    ggplot(phen.data,
           aes(Day_pi, Oocysts_feces, group=Mouse_ID, color=Mouse_strain)) +
    facet_wrap(~Infection_No)+
    geom_line() 


oocyst.summary <- ddply(phen.data, c("Mouse_strain", "Day_pi", "Infection_No"), summarize,
                      N    =  sum(!is.na(Oocysts_feces)),
                      Cmean = mean(Oocysts_feces, na.rm=TRUE),
                      Csd =   sd(Oocysts_feces, na.rm=TRUE),
                      Cse =   Csd/sqrt(N))


pdf("figures/Figure1a_oocystCounts.pdf", width=12, height=8)
all.oocysts.line +
    geom_line(data=oocyst.summary, 
              aes(Day_pi, Cmean, group=Mouse_strain, color=Mouse_strain), size=3)+
    geom_errorbar(data=oocyst.summary,
                  aes(x=Day_pi, y=Cmean, group=NULL,
                      ymin=Cmean-Cse, ymax=Cmean+Cse, color=Mouse_strain), size=2) +
    scale_x_continuous(limits=c(5, 12)) +
    theme_bw()
dev.off()

## qPCR data from Simone and Annica
## Ef18S qPCR data for different timepoints, 1st and 2nd infection and 3 biological replicates.
qpcr.df <- as.data.frame(read.csv("data/NMRI1st2ndqPCR_frSS_EfctRef_long.csv"))
long.qpcr <- melt(qpcr.df, measure.vars = c("rep1_ct", "rep2_ct", "rep3_ct"))
## Normalise to highest ct among all replicates
long.qpcr$norm.value <- (long.qpcr$value - max(long.qpcr$value))*-1


## Calculate sd and means on normalised values
stats.qpcr <- summaryBy(norm.value ~ dpi + inf,
          data = long.qpcr, 
          FUN = list(mean, min, max, sd, se))
stats.qpcr <- merge(stats.qpcr, long.qpcr)

################ STATUS ####################
## No transformation of data needed:
## plot on log2 scale: transformation on ggplot possible, but data is excluded 
## (due to neg values?): Negative values = more of target gene than of reference(?)

## Discussed Totta & Emanuel: EH will write now ~ 2 weeks
## he will focus more on upcoming grant application towards end of August
## Cool that we have so much qPCR data: find those genes in RNAseq data and 
## plot together

devSVG("~/Ef_RNAseq/figures/Figure1c_qPCR18S.svg")
qpcr18S <- ggplot(subset(stats.qpcr, stats.qpcr$gene %in% "Ef18S"), 
       aes(x = dpi, y = norm.value.mean, col = factor(inf))) +
  ggtitle("Eimeria 18S by qPCR") +
  labs(color = "Infection No.") +
  ## to get rid of label, rm aes(), but then size is not controlled. How fix this? ################ HERE ####################
  geom_errorbar(aes(ymin = norm.value.mean - norm.value.sd,
                    ymax = norm.value.mean + norm.value.sd,
                    width = 0.2)) +
  geom_point(aes(alpha = 0.7)) +
  geom_line() +
  ##  (delta-delta-ct): put in figure legend
    scale_y_continuous("Eimeria 18S fold change", 
                     labels = math_format(2^.x), 
                     breaks = c(0, 2, 4, 8, 16, 32),
                     limits = c(0, 32)) +
  scale_x_continuous("Day post infection", breaks = seq(0,8,1)) +
  theme_bw(20) +
  theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()




