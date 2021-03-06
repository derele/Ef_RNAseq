#### DIFFERENTIAL GENE EXPRESSION ANALYSIS
## needed from data import script are two raw count objects 
library(RSvgDevice)
library(statmod)
library(edgeR)
library(GGally)
library(reshape)
library(plyr)
library(gridExtra)
library(VennDiagram)
library(xtable)

## We know from 2_data_curation.R that we want to exclude:
excluded.samples <- c("NMRI_2ndInf_3dpi_rep1",
                      "NMRI_2ndInf_5dpi_rep2",
                      "NMRI_1stInf_0dpi_rep1")

All.RC <- as.matrix(read.table("output_data/RC_All_genes.csv", sep=","))

All.RC <- All.RC[, !colnames(All.RC)%in%excluded.samples]

## Mouse and Eimeria specific
Mm.RC <- All.RC[grepl('^ENSMUS.*', rownames(All.RC)),]
Ef.RC <- All.RC[grepl('^EfaB.*', rownames(All.RC)),]

pData <- read.table("output_data/Sample_pData.csv", sep=",")
pData <- pData[!pData$sample%in%excluded.samples, ]

## Mouse and Eimeria specific
Ef.pData <- pData[!pData$dpi%in%"0dpi", ]
Mm.pData <- pData[!pData$challenged%in%c("sporozoites", "oocysts"), ]

Ef.pData$grouped <- as.factor(as.character(Ef.pData$grouped))
Mm.pData$grouped <- as.factor(as.character(Mm.pData$grouped))

## Mouse and Eimeri spcific samples
Mm.RC <- Mm.RC[, colnames(Mm.RC)%in%Mm.pData$sample]
Ef.RC <- Ef.RC[, colnames(Ef.RC)%in%Ef.pData$sample] 

## Make sure the pData follows the same order as the sample count
## columns!!!
Mm.RC <- Mm.RC[, order(colnames(Mm.RC))]
Ef.RC <- Ef.RC[, order(colnames(Ef.RC))]

Mm.pData <- Mm.pData[order(Mm.pData$sample), ]
Ef.pData <- Ef.pData[order(Ef.pData$sample), ]



## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
Ef.design <- model.matrix(~0+Ef.pData$grouped)
colnames(Ef.design)  <- gsub("Ef.pData\\$grouped", "",
                             colnames(Ef.design))
## just to have the standard point delimitor
colnames(Ef.design)  <- gsub("_", ".",
                             colnames(Ef.design))

Ef.contrasts <- makeContrasts(
    ## just different stages in NMRI
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi, #1
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi,
    N3vsSpo = NMRI.1stInf.3dpi-NMRI.sporozoites,
    N3vsOoc = NMRI.1stInf.3dpi-NMRI.oocysts,
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi,
    N5vsSp = NMRI.1stInf.5dpi-NMRI.sporozoites,
    N5vsOo = NMRI.1stInf.5dpi-NMRI.oocysts,
    N7vsSp = NMRI.1stInf.7dpi-NMRI.sporozoites,
    N7vsOo = NMRI.1stInf.7dpi-NMRI.oocysts,
    SpvsOo=  NMRI.sporozoites-NMRI.oocysts,	#10
    ## differences in 1st vs 2nd infection
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi,	#11
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR5.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,	#15
    ## differences in mouse strains
    NvsB=NMRI.1stInf.5dpi-C57BL6.1stInf.5dpi,	#16
    NvsR=NMRI.1stInf.5dpi-Rag.1stInf.5dpi,
    BvsR=C57BL6.1stInf.5dpi-Rag.1stInf.5dpi,	#18
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),	#19
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi),	#21
##### the below is not working when Ef.1st.pass.model is produced yet ## intracellular versus extracellular stages 
    #NIntracell.vs.spo = ((NMRI.1stInf.3dpi - NMRI.1stInf.5dpi) + 
    #                      (NMRI.1stInf.3dpi - NMRI.1stInf.7dpi) + 
    #                       (NMRI.1stInf.5dpi - NMRI.1stInf.7dpi)) -
    #  (NMRI.sporozoites),
    #NIntracell.vs.oocy = () - (NMRI.oocysts),
    levels=Ef.design)


## use the alternate design like recommended by the EdgR manual (or
## eg. here: https://support.bioconductor.org/p/66952/) and specifiy
## interactions as contrasts!
Mm.design <- model.matrix(~0+Mm.pData$grouped)

colnames(Mm.design)  <- gsub("Mm.pData\\$grouped", "",
                             colnames(Mm.design))
## just to have the standard point delimitor
colnames(Mm.design)  <- gsub("_", ".",
                             colnames(Mm.design))

Mm.contrasts <- makeContrasts(
    ## just infections in NMRI !!! As we not have uninfected control
    ## using 2nd infected control here
    N3vs0 = NMRI.1stInf.3dpi-NMRI.2ndInf.0dpi, #1
    N5vs0 = NMRI.1stInf.5dpi-NMRI.2ndInf.0dpi, #2
    N7vs0 = NMRI.1stInf.7dpi-NMRI.2ndInf.0dpi, #3
    ## just infections in Black and Rag
    B5vs0 = C57BL6.1stInf.5dpi- C57BL6.1stInf.0dpi , #4
    R5vs0 = Rag.1stInf.5dpi - Rag.1stInf.0dpi,		#5
    ## differences in 1st vs 2nd infection
    N3.1stvsN3.2nd = NMRI.1stInf.3dpi-NMRI.2ndInf.3dpi, #6
    N5.1stvsN5.2nd = NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi,
    N7.1stvsN7.2nd = NMRI.1stInf.7dpi-NMRI.2ndInf.7dpi,
    B5.1stvsB5.2nd = C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi,
    R5.1stvsR5.2nd = Rag.1stInf.5dpi-Rag.2ndInf.5dpi,	#10
    ## differences in mouse strains
        RvsB=Rag.1stInf.0dpi-C57BL6.1stInf.0dpi,
	#NvsB=NMRI.0dpi-C57BL6.0dpi,
    ## differences in infection over time
    N3vsN5 = NMRI.1stInf.3dpi-NMRI.1stInf.5dpi,
    N5vsN7 = NMRI.1stInf.5dpi-NMRI.1stInf.7dpi,
    N3vsN7 = NMRI.1stInf.3dpi-NMRI.1stInf.7dpi,		#14
    ## what about differences in the reaction to first infection in
    ## different mouse strains??
    ##
    ## differences in 1st vs. 2nd depending on mouse strain
    N5.1stvsN52ndVSB5.1stvsB52nd = (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),		#15
    R5.1stvsR52ndVSB5.1stvsB52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (C57BL6.1stInf.5dpi-C57BL6.2ndInf.5dpi),
    R5.1stvsR52ndVSN5.1stvsN52nd = (Rag.1stInf.5dpi-Rag.2ndInf.5dpi)-
        (NMRI.1stInf.5dpi-NMRI.2ndInf.5dpi), 		#16
    levels=Mm.design)


get.my.models <- function (RC, cutoff, group,
                           design, contrasts,
                           norm="upperquartile"){
    keep <- rowSums(RC)>cutoff
    DGEList <- DGEList(RC[keep,], group=group)
    if(!is.null(norm)){
        DGEList <- calcNormFactors(DGEList, method=norm )
    }
    ## robust estimator https://support.bioconductor.org/p/65558/
    DGEList <- estimateDisp(DGEList, design = design, robust = TRUE)
    fit <- glmFit(DGEList, design) 
    fitLRT <- glmLRT(fit, contrast = contrasts)
    ALL.top <- topTags(fitLRT, 1000000)
    fitLRT.list <- lapply(colnames(contrasts),
                          function (x) glmLRT(fit,
                                              contrast = contrasts[,x]))
    names(fitLRT.list) <- colnames(contrasts)
    top.list <- lapply(fitLRT.list, function (x){
                           topTags(x, 1000000)[[1]]
                       })
    names(top.list) <- colnames(contrasts)
    gene.list <- lapply(top.list, function(x) {
                            rownames(x[x$FDR<0.01,]) ######### CHANGE FDR here
                        })
    return(list(ALL.top, top.list, gene.list, DGEList))
}

## filter: with 3000 as cutoff, bimodal distr. almost not visible (see
## A_data_curation.R)
Mm.1st.pass.model <- get.my.models(Mm.RC, cutoff=1000,
                                   group = Mm.pData$grouped,
                                   contrasts=Mm.contrasts,
                                   design=Mm.design,
                                   norm="upperquartile")

Ef.1st.pass.model <- get.my.models(Ef.RC, cutoff=10,
                                   group = Ef.pData$grouped,
                                   contrasts=Ef.contrasts,
                                   design=Ef.design,
                                   norm="upperquartile")


Mm.DE.test <- do.call(rbind, Mm.1st.pass.model[[2]])
Mm.DE.test$contrast <- gsub("(.*)\\.(ENSMUS.*?)$", "\\1", rownames(Mm.DE.test))
Mm.DE.test$gene <- gsub("(.*)\\.(ENSMUS.*?)$", "\\2", rownames(Mm.DE.test))
Mm.DE.test <- Mm.DE.test[, c("contrast", "gene", "logFC", "logCPM", "LR", "PValue", "FDR")]

Ef.DE.test <- do.call(rbind, Ef.1st.pass.model[[2]])
Ef.DE.test$contrast <- gsub("(.*)\\.(EfaB_.*?)$", "\\1", rownames(Ef.DE.test))
Ef.DE.test$gene <- gsub("(.*)\\.(EfaB_.*?)$", "\\2", rownames(Ef.DE.test))
Ef.DE.test <- Ef.DE.test[, c("contrast", "gene", "logFC", "logCPM", "LR", "PValue", "FDR")]


write.table(Mm.DE.test, "output_data/Mm_DEtest.csv", sep=",")
write.table(Ef.DE.test, "output_data/Ef_DEtest.csv", sep=",")


write.table(cpm(Mm.1st.pass.model[[4]]), "output_data/Mm_norm_counts.csv", sep=",")
write.table(cpm(Ef.1st.pass.model[[4]]), "output_data/Ef_norm_counts.csv", sep=",")


## how many genes are regulated:
melt(lapply(Mm.1st.pass.model[[3]], length))

melt(lapply(Ef.1st.pass.model[[3]], length))

## test for consistency with suppl file data:
by(Mm.DE.test, Mm.DE.test$contrast, function (x) table(x$FDR<0.01))
by(Ef.DE.test, Ef.DE.test$contrast, function (x) table(x$FDR<0.01))

Mm.FDR0.01 <- by(Mm.DE.test, Mm.DE.test$contrast, function (x) nrow(x[x$FDR<0.01, ]))
Mm.FDR0.01 <- cbind(Mm.FDR0.01)

Ef.FDR0.01 <- by(Ef.DE.test, Ef.DE.test$contrast, function (x) nrow(x[x$FDR<0.01, ]))
Ef.FDR0.01 <- cbind(Ef.FDR0.01)

test.FDR0.01 <- merge(Ef.FDR0.01, Mm.FDR0.01, by=0, all=TRUE)
rownames(test.FDR0.01) <- test.FDR0.01$Row.names
test.FDR0.01$Row.names <- NULL

test.FDR0.01 <- test.FDR0.01[rowSums(test.FDR0.01, na.rm=TRUE)>10, ] 

test.FDR0.01 <- test.FDR0.01[order(test.FDR0.01$Mm.FDR0.01, test.FDR0.01$Ef.FDR0.01,
                               decreasing=TRUE), ]

test.FDR.tabl <- xtable(test.FDR0.01, digits=c(NA, 0, 0))

print(test.FDR.tabl,
      type = "html", file = "tables/Table2_DE_tests.html", include.rownames = T,
      format.args = list(big.mark = ",", decimal.mark = "."))

########### MOUSE #################
Mm.infection.difference <-
    venn.diagram(Mm.1st.pass.model[[3]][c("N3vs0", "N5vs0", "N7vs0")],
                 filename=NULL)

devSVG("figures/Figure2a_vennMmInfection.svg")
grid.draw(Mm.infection.difference)
dev.off()

########## EIMERIA ##################

Ef.infection.LCInteral.difference <-
    venn.diagram(Ef.1st.pass.model[[3]][c("N3vsN5", "N3vsN7", "N5vsN7")],
                 filename=NULL)

devSVG("figures/Figure3a_vennEfCycleInternal.svg")
grid.draw(Ef.infection.LCInteral.difference)
dev.off()

###
groups <- apply(Mm.contrasts[, 1:14], 2, function (x) names(x[x!=0]))
groups <- gsub("\\.", "_", groups)

pdf("Supplement/allMA.pdf", width=12, height=18)
par(mfrow=c(7,4))
sapply(colnames(groups), function (x){
    plotSmear(Mm.1st.pass.model[[4]],
              pair = c(groups[1, x], groups[2, x]), de.tags=Mm.1st.pass.model[[3]][[x]],
              main ="Mouse")})
groups <- apply(Ef.contrasts[, 1:14], 2, function (x) names(x[x!=0]))
groups <- gsub("\\.", "_", groups)
sapply(colnames(groups), function (x){
    plotSmear(Ef.1st.pass.model[[4]],
              pair = c(groups[1, x], groups[2, x]), de.tags=Ef.1st.pass.model[[3]][[x]],
              main = "E. falciformis")})
dev.off()
