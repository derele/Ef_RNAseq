## DIFFERENTIAL GENE EXPRESSION ANALYSIS ON EIMERIA SEQS #####

#setwd("/home/totta/Ef_RNAseq")

if(!exists("raw.counts.4.bg")){
	source("2_edgeR_diff.R")
        }

Ef.RC <- raw.counts.4.bg(Ef.bg)

pdf("figures/rep_pairs_EfG.pdf", width = 70, height=70) # So far not working - because R interrupt?
ggpairs(Ef.RC[[3]])
dev.off()

## use the alternate designE like .... the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
#MOUSE: design <- model.matrix(~0+pData(Ef.bg)$grouped)

##################### For EIMERIA: #################################

## Summary table of transcript counts per sample
## MOUSE table
tableM <- as.data.frame(colSums(mouse.RC[[3]]))
#colnames(tableM)[,2] <- "Read.number"
## MOUSE w/o day 0 
tableM$Samples <- rownames(tableM)
library(stringr)
tableMsmall <- tableM[!(str_detect(rownames(tableM),"^.*0dpi.*$")), ]

## EIMERIA table
tableEf <- as.data.frame(colSums(Ef.RC[[3]]))
#colnames(tableEf)[,1] <- "Read.number" - doesnt work?
#tableEf$Samples <- colnames(Ef.RC[[3]])
#print.xtable(tableEf, type = "latex", file = "tableEf-reads.tex")
tableEf$Samples <- rownames(tableEf)
## EIMERIA table w/o oocysts and sporozoites, for mouse comparison
tableEfsmall <- tableEf[!(str_detect(tableEf$Samples, "^NMRI_(oocysts|sporozoites)$")), ]

## ADD EIMERIA % to EIMERIAsmall table
tableEfsmall$percent.Ef.reads <- (tableEfsmall[,1]/tableMsmall[,1])*100

library(scales)
ggplot()+
	geom_point(data = tableEf, aes(x = Samples, y = as.numeric(colSums(Ef.RC[[3]])))) +
	scale_y_log10(labels = comma) +
	ylab("Number of reads") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

#geom_point(data = tableM, aes(x = , y = as.numeric(colSums(Ef.RC[[3]])))

##############################
## day 0 rm from pData and designE
designE <- model.matrix(~0+pData(Ef.bg)$grouped)# 0+ in the model formula is an
						# instruction not to include an intercept column and
						# instead to include a column for each group

colnames(designE)  <- gsub("pData\\(Ef.bg\\)\\$grouped", "",
                             colnames(designE))
## just to have the standard point delimitor
colnames(designE)  <- gsub("_", ".",
                             colnames(designE))
#### lower number for Ef, higher for mouse?? Play
keep <- rowSums(Ef.RC[[3]])>5 ## filter lowly expressed raw counts
GM.E <- DGEList(Ef.RC[[3]][keep,], group = pData(Ef.bg)$grouped)

GM.E <- calcNormFactors(GM.E)

## robust estimator might help
## https://support.bioconductor.org/p/65558/
GM.E <- estimateDisp(GM.E, design = designE, robust = TRUE)

pdf("figures/rep_pairs_NORM_EfG.pdf", width = 70, height=70)
ggpairs(cpm(GM.E)) #ggpairs makes matrix of plots
dev.off()

pdf("figures/Ef_genes_mds.pdf")
plotMDS(GM.E, )
dev.off()

pdf("figures/Ef_genes_hclust.pdf")
plot(hclust(dist(t(cpm(GM.E)))))
dev.off()

fitE <- glmFit(GM.E, designE) # Is it correct to take the dispersions as input here??

my.contrastsE <- makeContrasts(
    ## differences in 1st vs 2nd infection, (0dpi excluded)
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi, # 10 returned from topTags
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi, # nothing
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi, # nothing
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi, # nothing
    R5.1stvsR2.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi, # nothing
    ## differences in infection over time
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi, # 47 returned
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi, # 2684 returned
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi, # 2227 returned -> day 3 and 7 more alike than 5 and 7?
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi), # nothing
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),  # nothing. 
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi), # nothing
    levels=designE)

fitLRT <- glmLRT(fitE, contrast = my.contrastsE) # when contrasts are given,
						 # tests the likelihood that that two contrasts are equal (H-null)

ALL.top <- topTags(fitLRT, 1000000)

fitLRT.list <- lapply(colnames(my.contrastsE),
                          function (x) glmLRT(fitE, contrast = my.contrastsE[,x]))
names(fitLRT.list) <- colnames(my.contrastsE)

top.listE <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
names(top.listE) <- colnames(my.contrastsE)

######## How to deal with duplicates: _1, _2, _3..? Searched for Tophat,
	#  Bowtie and Cufflinks add suffix to id(gene_id/_id but did not find explanation yet.
gene.list.E <- lapply(top.listE, function(x) {
                            rownames(x[x$FDR<0.05,])
                        })

