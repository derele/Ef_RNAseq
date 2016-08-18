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
library(RColorBrewer)

setwd("~/Ef_RNAseq/")
phen.data <- read.csv("data/Oocysts_output_weight_SS_longdata.csv")

## here just a hack to take the first occurence of a ID at a time point
phen.data <- phen.data[!duplicated(phen.data[, c("Day_pi", "Mouse_ID", "Infection_No")]), ]

#all.oocysts.line <-
#    ggplot(phen.data,
#           aes(Day_pi, Oocysts_feces, group=Mouse_ID, color=Mouse_strain)) +
#    facet_wrap(~Infection_No)+
#    geom_line() 

all.oocysts.line <-
  ggplot(phen.data,
         aes(Day_pi, Oocysts_feces, color=Mouse_strain)) +
  facet_wrap(~Infection_No)

  #geom_line() 


oocyst.summary <- ddply(phen.data, c("Mouse_strain", "Day_pi", "Infection_No"), summarize,
                      N    =  sum(!is.na(Oocysts_feces)),
                      Cmean = mean(Oocysts_feces, na.rm=TRUE),
                      Csd =   sd(Oocysts_feces, na.rm=TRUE),
                      Cse =   Csd/sqrt(N))

## create colors for plot
palette.colors <- brewer.pal(8, "Dark2")
my.colors.ooc <- palette.colors[c(1,2,8)]

pdf("figures/Figure1a_oocystCounts.pdf", width=12, height=8)
all.oocysts.line +
    geom_point(data=oocyst.summary, 
               aes(Day_pi, Cmean)) +
    geom_line(data=oocyst.summary, 
              aes(Day_pi, Cmean),
              size = 0.7) +
    geom_errorbar(data=oocyst.summary,
                  aes(x=Day_pi, y=Cmean, #group=NULL,
                      ymin=Cmean-Cse, 
                      ymax=Cmean+Cse), 
                  size=1,
                  width = 0.3) +
    scale_color_manual(values = my.colors.ooc) +
    scale_y_continuous("Oocyst number in feces", 
                     labels = comma, 
                     breaks = c(0, 1000000, 2000000, 3000000, 4000000)) +
    scale_x_continuous("Day post infection",
                       breaks = c(0, 3, 5, 7, 9, 11, 13, 15)) + 
    theme_bw(20) +
    theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    ggtitle("Oocyst counts in first and second infection")
## not saving the latter part - only part creating the plot   !!!!!!!!!!!!!!!!!!!
#ggsave(file = "figures/Figure1a_oocystsCounts.svg", plot = all.oocysts.line)
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
          FUN = list(mean, min, max, sd))
stats.qpcr <- merge(stats.qpcr, long.qpcr)

################ STATUS ####################
## Cool that we have so much qPCR data: find those genes in RNAseq data and 


## plot together
## create colors for plot
palette.colors2 <- brewer.pal(11, "BrBG")
my.colors.qpcr <- palette.colors2[c(3,9)]

##### plot goves error message due to my.colors code line, but plot looks fine...
#devSVG(file = "figures/Figure1c_qPCR18S.svg")
qpcr18S <- ggplot(subset(stats.qpcr, stats.qpcr$gene %in% "Ef18S"), 
       aes(x = dpi, y = norm.value.mean, col = factor(inf))) +
  ggtitle("Eimeria 18S by qPCR") +
  labs(color = "Infection No.") +
  ## to get rid of label, rm aes(), but then size is not controlled. How fix this? ################ HERE ####################
  geom_errorbar(aes(ymin = norm.value.mean - norm.value.sd,
                    ymax = norm.value.mean + norm.value.sd,
                    width = 0.3)) +
  geom_point(aes(alpha = 0.7)) +
  geom_line(size = 0.7) +
  scale_color_manual(values = my.colors.qpcr) +
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
ggsave(file = "figures/Figure1c_qPCR18S.svg", plot = qpcr18S)
dev.off()




