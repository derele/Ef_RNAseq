## Reading and plotting Simones data on oocyst output and weigh change upon Eimeria infection
## of NMRI, C57BL/6 and Rag1-/- mice
library(ggplot2)
library(gridExtra)

setwd("~/Ef_RNAseq/")
phen.data <- read.csv("data/Oocysts_output_weight_SS_longdata.csv")

## Immunocompetent mice oocyst output, 1st adn 2nd infection
bl6NMRI.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% c("NMRI", "C57BL6")),
       aes(x = Day_pi, y = Oocysts_feces, col = factor(Infection_No))) +
  geom_point(position = "jitter") +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Oocyst count in feces") +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
  scale_colour_discrete("Infection") +
  ggtitle("Oocysts output in immunocompetent mice")
## In immunodeficient mice, 1st and 2nd infection
rag.oocysts <- ggplot(subset(phen.data, phen.data$Mouse_strain %in% "Rag"),
       aes(x = Day_pi, y = Oocysts_feces, col = factor(Infection_No))) +
  geom_point(position = "jitter") +
  theme_bw(20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Oocyst count in feces") +
  scale_x_continuous("Day post infection", breaks = seq(0, 15, 1)) + 
  scale_colour_discrete("Infection") +
  ggtitle("Oocysts output in immunocompromised mice")
rag.oocysts  

ggsave("figuresANDmanuscript/immuneVSrag_oocysts.png")
grid.arrange(bl6NMRI.oocysts, rag.oocysts, ncol = 1)

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