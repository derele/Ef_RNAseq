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

## Reshape so that oocysts and weights are variables in same column
## (to enable calculation of means and sd later on)
long.phen.data <- melt(phen.data, measure.vars = c("Oocysts_feces", "Normalized_weight"))
## Calculate means and sd with dplyr
#phen.stats <- ddply(long.phen.data, .(Mouse_ID, 
#                                      variable,
#                                      value), function(x){
#                     c(mean = mean(x$value), sd = sd(x$value), CI = confint(x$value))
#                     })

## Create Infection number and mouse-group category to be able to identify replicates
long.phen.data$id <- paste(phen.data$Day_pi,
                           #phen.data$Groups_unique,
                           phen.data$Mouse_strain,
                           long.phen.data$variable,
                           phen.data$Infection_No,
                           sep = "")
phen.summary <- summaryBy(value ~ id, 
                          data = long.phen.data[], 
          FUN = list(mean, min, max, sd))
phen.summary$Day_pi <- as.numeric(str_match(phen.summary$id, "^(\\d*).*|\\d*.*\\d$")[,2])
phen.summary$Inf_number <- str_match(phen.summary$id, "^\\d*.*(\\d)$")[,2]
phen.summary$Mouse_strain <- str_match(phen.summary$id, "^\\d*(NMRI|Rag|C57BL6).*$")[,2]
phen.summary$Phenotype <- str_match(phen.summary$id, "^\\d*.*(Normalized_weight|Oocysts_feces)\\d$")[,2]

## Make function for plotting
## dat = dataframe, 
## strain = mice strain(s) to be used in plot
## inf = infection number, ie 1 or 2
## phen = phenotype to be plotted, ie. weight or oocysts
#drawplot <- function(dat, phen) {
#  ggplot(subset(dat, dat$Phenotype %in% phen),
#       aes(x = Day_pi, y = value.mean, col = factor(Mouse_strain))) +
#  geom_line(linetype = Infection_No) +
#  geom_point() +
#  theme_bw(20) +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), 
#        axis.line = element_line(colour = "black"),
#        axis.title.y = element_blank(),
#        legend.key = element_blank()) +
#  scale_y_continuous("Average weight as % n\ of pre-infection", 
#                     labels = comma) +#, 
#                     #breaks = seq(75, 130, 10),
#                     #limits = c(75, 130)) +
#  #scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
#  #scale_colour_discrete("Infection") +
#  ggtitle("Weight change by immune status")
#  }

#drawplot(phen = phen.summary, phen = "Oocysts_feces")

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

limits <- aes(ymax = phen.summary$value.mean + phen.summary$value.sd,
              ymin = phen.summary$value.mean - phen.summary$value.sd)


weight.averages <- ggplot(subset(phen.summary,phen.summary$Phenotype %in% "Normalized_weight"),
                          aes(x = Day_pi, y = value.mean, col = factor(Mouse_strain))) +
  geom_errorbar(x = Day_pi, ymax = phen.summary$value.mean + phen.summary$value.sd,
                    ymin = phen.summary$value.mean - phen.summary$value.sd) +
  #geom_line(aes(linetype = Inf_number)) +
  geom_point(aes(symbols = Inf_number)) +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Average weight as % n\ of pre-infection", 
                     labels = comma) + #, 
                     #breaks = seq(75, 130, 10),
                     #limits = c(75, 130)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
  scale_colour_discrete("Strain") +
  ggtitle("Weight change by immune status")

ggsave("figuresANDmanuscript/weight_averages.png", plot = weight.averages)

oocysts.averages <- ggplot(subset(phen.summary,
                                  phen.summary$Phenotype %in% "Oocysts_feces"),
                           aes(x = Day_pi, y = value.mean, col = factor(Mouse_strain))) +
  geom_errorbar(aes(ymin = value.min, ymax = value.max)) +
  geom_line(aes(linetype = Inf_number)) +
  geom_point() + 
#  facet_grid(Infection_No ~ Mouse_strain) +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Oocyst counts in feces", 
                     labels = comma) + #,
                     #breaks = seq(0, 6000000, 2000000),
                     #limits = c(0, 4000000)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 3)) + 
  scale_colour_discrete("Strain") +
  ggtitle("Parasite load by immune status")

ggsave("figuresANDmanuscript/oocyst_averages.png", plot = oocysts.averages)


## Weight change in immune competent and compromised (respectively) in panel
bl6NMRI_1.weight <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% c("NMRI", "C57BL6")
                                & phen.data$Infection_No %in% 1),
      aes(x = Day_pi, y = Normalized_weight, col = factor(Mouse_ID))) +
  geom_line() +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Weight in % compared to pre infection", 
                     labels = comma, 
                     breaks = seq(75, 130, 10),
                     limits = c(75, 130)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
  scale_colour_discrete("Infection") +
  ggtitle("Weight change in immunocompetent mice")
## Immuno compromised Rag mice weight change
rag.weight <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% "Rag"),
                         aes(x = Day_pi, y = Normalized_weight, col = factor(Mouse_ID))) +
  geom_line() +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Weight in % compared to pre infection", 
                     labels = comma, 
                     breaks = seq(75, 130, 10),
                     limits = c(75, 130)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
  scale_colour_discrete("Infection") +
  ggtitle("Weight change in immunocompromised mice")

immuneVSrag.weight <- grid.arrange(bl6NMRI.weight, rag.weight, 
                                   ncol = 1,
                                   left = textGrob("Weight in % compared to control", 
                                                   rot = 90, 
                                                   vjust = 1,
                                                   gp = gpar(fontsize = 20)))
ggsave("figuresANDmanuscript/immuneVSrag_weight.png", plot = immuneVSrag.weight)


## Immunocompetent mice oocyst output, 1st adn 2nd infection
bl6NMRI_1.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% c("NMRI", "C57BL6") 
                                 & phen.data$Infection_No %in% 1),
       aes(x = Day_pi, y = Oocysts_feces, group = Mouse_ID, col = factor(Mouse_strain))) +
  geom_point() +
  theme_bw(16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Oocyst count in feces", 
                     labels = comma, 
                     breaks = seq(0, 6000000, 1000000),
                     limits = c(0, 4000000)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 3), limits = c(0, 15)) + 
  scale_colour_discrete("Infection", guide_legend(F)) +
  ggtitle("Competent, 1st")

## 2nd infection
bl6NMRI_2.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% c("NMRI", "C57BL6") 
                                 & phen.data$Infection_No %in% 2),
                          aes(x = Day_pi, y = Oocysts_feces, group = Mouse_ID, col = factor(Mouse_strain))) +
  geom_point() +
  theme_bw(16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous("Oocyst count in feces", labels = comma, 
                     breaks = seq(0, 6000000, 1000000), 
                     limits = c(0, 4000000)) +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 3), limits = c(0, 15)) + 
  scale_colour_discrete("Infection", guide_legend(F)) +
  ggtitle("Competent, challenge")

## Combine into 1 plot
oocystImmune_panel <- grid.arrange(bl6NMRI_1.oocysts, bl6NMRI_2.oocysts, 
             ncol = 2)

##################################################
## In immunodeficient mice, 1st and 2nd infection
rag_1.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% "Rag"
                               & phen.data$Infection_No %in% 1),
       aes(x = Day_pi, y = Oocysts_feces, group = Mouse_ID, col = factor(Mouse_strain))) +
  geom_line() +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 3), limits = c(0, 15)) + 
  scale_y_continuous(labels = comma, 
                     breaks = seq(0, 6000000, 1000000), 
                     limits = c(0, 4000000)) +
  theme_bw(16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank()) +
  scale_colour_discrete("Infection", guide_legend(F)) +
  ggtitle("Compromised, 1st")

rag_2.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% "Rag"
                               & phen.data$Infection_No %in% 2),
                        aes(x = Day_pi, y = Oocysts_feces, group = Mouse_ID, col = factor(Mouse_strain))) +
  geom_line() +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 3), limits = c(0, 15)) + 
  scale_y_continuous(labels = comma, breaks = seq(0, 6000000, 1000000), limits = c(0, 4000000)) +
  theme_bw(16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank()) +
  scale_colour_discrete("Infection", guide_legend(F)) +
  ggtitle("Compromised, challenge")

oocystsRag_panel <- grid.arrange(rag_1.oocysts, rag_2.oocysts,
             ncol = 2)

immuneVSrag.oocysts <- grid.arrange(oocystImmune_panel, oocystsRag_panel, 
                                    ncol = 1,
                                    top = textGrob("Oocyst counts in different mice",
                                                   vjust = 0.4,
                                                   hjust = 0.5,
                                                   gp = gpar(fontsize = 25)),
                                    bottom = textGrob("Day post infection",
                                                      vjust = 0.4,
                                                      hjust = 0.3,
                                                      gp = gpar(fontsize = 25)),
                                    left = textGrob("Oocyst counts", 
                                                    rot = 90, 
                                                    vjust = 1,
                                                    gp = gpar(fontsize = 25)))
ggsave("figuresANDmanuscript/immuneVSrag_oocysts.png", plot = immuneVSrag.oocysts)

## Plot weight - 1st infection - as function of time post infection
weight.plot <- ggplot(subset(phen.data, as.numeric(phen.data$Infection_No) %in% 1), 
                      aes(x = Day_pi, y = Normalized_weight, col = factor(Mouse_strain))) +
  geom_point() +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) +
  ggtitle("Weight change in first infection") +
  theme_bw(20)
weight.plot
dev.off()

## Plot oocysts in feces - 1st infection - as function of time post infection
oocyst.plot <- ggplot(subset(phen.data, as.numeric(phen.data$Infection_No) %in% 1), 
                      aes(x = Day_pi, y = Oocysts_feces, col = factor(Mouse_strain))) +
  #geom_point(position = "jitter") +
  geom_point() +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) +
  ggtitle("Oocyst output over course of infection") +
  theme_bw(20)
#  scale_y_log10()
oocyst.plot
dev.off()

## Nothing really from this comparison:
oocystVSweight <- ggplot(subset(phen.data, as.numeric(phen.data$Infection_No) %in% 1), 
                      aes(x = Oocysts_feces, y = Normalized_weight, col = factor(Mouse_strain))) +
  #geom_point(position = "jitter") +
  geom_point()  +
  #scale_x_continuous("Oocysts in feces")#, breaks = seq(0, 15, 1))
  scale_x_log10() + 
  theme_bw(20)
oocystVSweight
dev.off()