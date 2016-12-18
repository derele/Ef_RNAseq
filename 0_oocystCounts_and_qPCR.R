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


## overall output for statistical tests
## reduce to 10dpi to have all strains up to the same day
phen.data[phen]

oocyst.sums <- ddply(phen.data, c("Mouse_strain", "Mouse_ID", "Infection_No"), summarize,
                     N    =  sum(!is.na(Oocysts_feces)),
                     Csum    =  sum(Oocysts_feces))

wilcox.test(oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"C57BL6"&
                             oocyst.sums$Infection_No%in%"1"],
            oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"Rag"&
                             oocyst.sums$Infection_No%in%"1"]
            )    
## ns

wilcox.test(oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"C57BL6"&
                             oocyst.sums$Infection_No%in%"2"],
            oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"Rag"&
                             oocyst.sums$Infection_No%in%"2"]
            )    
## ns
wilcox.test(oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"C57BL6"&
                             oocyst.sums$Infection_No%in%"1"],
            oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"C57BL6"&
                             oocyst.sums$Infection_No%in%"2"]
            )    
## significant !!!

wilcox.test(oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"Rag"&
                             oocyst.sums$Infection_No%in%"1"],
            oocyst.sums$Csum[oocyst.sums$Mouse_strain%in%"Rag"&
                             oocyst.sums$Infection_No%in%"2"]
            )    
## Rag has also a reduced output in secondary !!! not significant
## though... but you know you can't statistically show absence of diff



oocyst.summary <- ddply(phen.data, c("Mouse_strain", "Day_pi", "Infection_No"), summarize,
                        N    =  sum(!is.na(Oocysts_feces)),
                        Cmean = mean(Oocysts_feces, na.rm=TRUE),
                        Csd =   sd(Oocysts_feces, na.rm=TRUE),
                        Cse =   Csd/sqrt(N))

## create colors for Figure 1a
palette.colors <- brewer.pal(8, "Dark2")
my.colors.ooc <- palette.colors[c(1,2,8)]
#pdf("figures/Figure1a_oocystCounts.pdf", width=12, height=8)

#labels for facet grid for Figure 1a
my.labels <- c("1" = "1st infection", "2" = "2nd infection")

## Plotting Figure 1a
all.oocysts.line <- ggplot(phen.data,
         aes(Day_pi, Oocysts_feces, color=Mouse_strain)) +
  facet_wrap(~Infection_No, labeller = labeller(Infection_No = my.labels)) + 
  geom_point(data=oocyst.summary, aes(Day_pi, Cmean)) +
  geom_line(data=oocyst.summary, aes(Day_pi, Cmean), size = 0.7) +
  geom_errorbar(data=oocyst.summary, aes(x=Day_pi, y=Cmean, ymin=Cmean-Cse, ymax=Cmean+Cse), 
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
      panel.grid.minor = element_blank())
  #ggtitle("Oocyst counts in first and second infection")
ggsave(file = "figures/Figure1a_oocystsCounts.svg", height = 8, width = 12, plot = all.oocysts.line)
#dev.off()

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
## plot together?

## create colors for Figure 1b
palette.colors2 <- brewer.pal(11, "BrBG")
my.colors.qpcr <- palette.colors2[c(3,9)]
my.ylab <- expression(paste(italic("Eimeria"), " 18S fold change by qPCR"))

##### Plotting figure 1b
qpcr18S <- ggplot(subset(stats.qpcr, stats.qpcr$gene %in% "Ef18S"), 
       aes(x = dpi, y = norm.value.mean, col = factor(inf))) +
  #ggtitle("Eimeria 18S by qPCR") +
  labs(color = "Infection No.") +
  ## to get rid of label, rm aes(), but then size is not controlled. How fix this? ################ HERE ####################
  geom_errorbar(aes(ymin = norm.value.mean - norm.value.sd,
                    ymax = norm.value.mean + norm.value.sd,
                    width = 0.3)) +
  geom_point() + #aes(alpha = 0.7)) +
  geom_line(size = 1) +
  scale_color_manual(values = my.colors.qpcr) +
  ##  (delta-delta-ct): put in figure legend
  scale_y_continuous(my.ylab, 
                     labels = math_format(2^.x), 
                     breaks = c(0, 2, 4, 8, 16, 32),
                     limits = c(0, 32)) +
  scale_x_continuous("Day post infection", breaks = seq(0,8,1)) +
  theme_bw(32) + #needs to be bigger than oocyst plot
  theme(legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file = "figures/Figure1b_qPCR18S.svg", height = 12, width = 11, plot = qpcr18S)
#dev.off()




