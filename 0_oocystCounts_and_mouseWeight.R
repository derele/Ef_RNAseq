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


## Do we really need the weight data? If yes do as above.
