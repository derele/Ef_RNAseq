### Check distributions of transcript data

###########################
library(reshape2)
############## MOUSE ##################
## rearrange GM.E for ggplot format
long.Mcounts <- melt(GM$counts)
long.counts$strain <- gsub("_.*", "", long.Mcounts$Var2)
long.Mcounts$strain <- gsub("_.*", "", 
		long.Mcounts$Var2)

long.normMcounts <- melt(cpm(GM))
long.normMcounts$strain <- gsub("_.*", "", 
		long.normMcounts$Var2)

## RAW MOUSE DISTRIBUTION
pdf("figures/distributions_mouse5000.pdf", 
	width = 12, height = 7)
ggplot(long.Mcounts, 
	aes(value)) + 
	stat_density(geom="line") +
	scale_x_log10("Read counts (log10)",
			labels = trans_format("log10",
			math_format(10^.x))) +
	facet_wrap(~Var2) +
	ggtitle("Mouse, non-norm., cutoff = 5000 reads") +
	ylim(0, 0.6)
dev.off()

## NORMALISED MOUSE DISTRIBUTION
pdf("figures/distributions_normmouse5000.pdf", 
	width = 12, height = 7)
ggplot(long.normMcounts, 
	aes(value)) + 
	stat_density(geom="line") +
	scale_x_log10("Read counts (log10)",
			labels = trans_format("log10",
			math_format(10^.x))) +
	facet_wrap(~Var2) +
	ggtitle("Mouse, norm., cutoff = 5000 reads") +
	ylim(0, 0.6)
dev.off()

############## PARASITE ###################
long.Ecounts <- melt(GM.E$counts)
long.Ecounts$strain <- gsub("_.*", "", long.Ecounts$Var2)
long.Ecounts$strain <- gsub("_.*", "", 
		long.Ecounts$Var2)

long.normEcounts <- melt(cpm(GM.E))
long.normEcounts$strain <- gsub("_.*", "", 
		long.normEcounts$Var2)

excludedE <- long.normEcounts[long.normEcounts$Var2%in%
		c("NMRI_2ndInf_7dpi_rep2",
		"NMRI_2ndInf_5dpi_rep2",
		"NMRI_2ndInf_3dpi_rep1",
		"NMRI_2ndInf_3dpi_rep2"),]
## See counts for these samples; there are very few Ef trascr. colSums(GM.E$counts)[order(colSums(GM.E$counts))]OR
# table(GM.E$counts[,"NMRI_2ndInf_3dpi_rep1" ]>20)summary(GM.E$counts[,"NMRI_2ndInf_3dpi_rep1" ])

## look at linear x-scale
## RAW EF DATA
pdf("figures/distributions_Ef5.pdf", 
	width = 12, height = 7)
ggplot(long.Ecounts, 
	aes(value)) +
		stat_density(geom="line") +
		scale_x_log10("Read counts (log10)", 
				labels = trans_format("log10", 
				math_format(10^.x))) +
		facet_wrap(~Var2) +
		ggtitle("Eimeria, non-norm., cutoff = 5 reads")
dev.off()

## NORMALISED EF DATA DISTRIBUTION
pdf("figures/distributions_Efnorm5.pdf", 
	width = 12, height = 7)
ggplot(long.normEcounts, 
	aes(x = value))+ 
		stat_density(geom="line")+
		scale_x_log10("Read counts (log10)",
				labels = trans_format("log10", 
				math_format(10^.x))) +
		facet_wrap(~Var2)+
		ggtitle("Eimeria, norm., cutoff = 5 reads") +
		ylim(0, 5)
dev.off()

## Plot of excluded distributions only
pdf("figures/distributions_Efexcluded.pdf", 
	width = 10, height = 7)
ggplot(excludedE, 
	aes(value, ..density.. , 
	color = Var2)) +
	geom_density() +
	scale_x_log10() +
	facet_wrap(~Var2)
dev.off()

#pdf("figures/distributions_Ef50000.pdf", 
#	width = 12, height = 7)
#ggplot(long.Ecounts, 
#	aes(value, ..density.. , 
#		color = Var2)) +
#		geom_density() +
#		scale_x_log10() +
#		facet_wrap(~Var2)
#dev.off()
### original plot
#pdf("figures/distributions_Efnorm50000.pdf", 
#	width = 12, height = 7)
#ggplot(long.normEcounts, 
#	aes(value, ..density..))+ 
#		#color = Var2)) +
#		guides(fill=F) +
#		geom_density() +
#		scale_x_log10("Read counts (log10)") +
#		facet_wrap(~Var2)
#dev.off()

######## COMMENTS ##########
# We should probably remove 2 or 4 samples in "excludedE" object
# and re-run analysis and DEG without them.
# These are the four samples with least transcript numbers, which
# probably explains the weird distributions. 
