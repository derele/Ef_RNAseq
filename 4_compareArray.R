## RNAseq - Microarray Comparison
## read the array data as done for Schmidt et al.  "Eimeria
## falciformis infection of the mouse caecum identifies opposing roles
## of IFNg-regulated host pathways for the parasite development"
## DOI:10.1038/mi.2013.115

#source("https://bioconductor.org/biocLite.R")
# Check checkpoint package for package control.
#if(!require(Biobase)) biocLite("Biobase") # Imports package if user does not have it
#if(!require(RSvgDevice)) biocLite("RSvgDevice") # Imports package if user does not have it
#if(!require(GGally)) biocLite("GGally") # Imports package if user does not have it
#if(!require(ggplot2)) biocLite("ggplot2") # Imports package if user does not have it
## For below package does not always work: then manual import is necessary
#if(!require(mgug4122a.db)) biocLite("mgug4122a.db") # Imports package if user does not have it

library(mgug4122a.db)
library(Agi4x44PreProcess) # !!!! NOT available for R 3.3.0 - NOT on Bioconductor anymore
library(Biobase)
library(RSvgDevice)
library(GGally)
library(ggplot2)

if(!exists("Mm.1st.pass.model")){
    source("2_edgeR_diff.R")
}

if(!exists("gene2GO")){
    source("3_annotations.R")
}

#targets <- readTargets("/data/Mouse_arrays/2012_arrays/targets.txt", sep=";")
targets <- readTargets("/SAN/Mouse_arrays/2012_arrays/targets.txt", sep=";")

RG <- read.maimages(targets, path= "/SAN/Mouse_arrays/2012_arrays/", source = "agilent", 
                    other.columns = list(IsFound = "gIsFound",
                      IsWellAboveBG = "gIsWellAboveBG", 
                      IsSaturated = "gIsSaturated",
                      IsFeatNonUnifOF = "gIsFeatNonUnifOL", 
                      IsFeatPopnOL = "gIsFeatPopnOL",
                      ChrCoord = "chr_coord"), 
                    columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal",
                      R = "rMedianSignal",
                      Rb = "rBGMedianSignal"), verbose = T, 
                    sep = "\t", quote = "")

array.names <- apply(targets, 1, function (x) paste(x[c(2, 4,5)], collapse="."))
colnames(RG) <- array.names
rownames(targets) <- array.names

RG <- backgroundCorrect(RG, method="minimum", offset=1)
MA <- normalizeWithinArrays(RG, method="loess", iterations=10)
MA <- normalizeBetweenArrays(MA, method="Aquantile")
                             
RG.n <- RG.MA(MA)
RG.n$other <- RG$other

RG.f <- filter.probes(RG.n,
                      control=TRUE,
                      wellaboveBG=TRUE,
                      isfound=TRUE,
                      wellaboveNEG=TRUE,
                      sat=TRUE,
                      PopnOL=TRUE,
                      NonUnifOL=T,
                      nas=TRUE,
                      limWellAbove=75,
                      limISF=75,
                      limNEG=75,
                      limSAT=75,
                      limPopnOL=75,
                      limNonUnifOL=75,
                      limNAS=100,
                      makePLOT=F,annotation.package="mgug4122a.db",
                      flag.counts=T, targets)

RG.p <- summarize.probe(RG.f, makePLOT=FALSE, targets)
MA.p <- MA.RG(RG.p)

ES.p <- build.eset(RG.p, targets, makePLOT=FALSE,
                   annotation.package="mgug4122a.db")
## Define the dye-swap as biological replicates!!!
biolrep <- targets$SampleNumber

design <- modelMatrix(targets, ref="Uninf")
cont.matrix <- makeContrasts("X24h" , "X144h", "X24h-X144h",
                             levels=design)

MA.g <- avereps(MA.p, ID=MA.p$genes$GeneName)
RG.g <- RG.MA(MA.g)
ES.g <- build.eset(RG.g, targets, makePLOT=FALSE,
                   annotation.package="mgug4122a.db")

corfit.g <- duplicateCorrelation(MA.g, design=design, ndups=1, block=biolrep)

fit.g <- lmFit(MA.g, design=design, block=biolrep, cor=corfit.g$consensus)

fit2.g <-  contrasts.fit(fit.g, cont.matrix)
fit2.g <- eBayes(fit2.g)

Array.logFC <- topTable(fit2.g, n=400000)[, c("X24h", "X144h")]

RNAseq.logFC <- Mm.1st.pass.model[[1]][, c("logFC.N3vs0", "logFC.N5vs0",
                                           "logFC.N7vs0", "logFC.B5vs0",
                                           "logFC.R5vs0" )]

RNAseq.logFC <- merge(RNAseq.logFC, annot.frame, by.x = 0, by.y = "ensembl_id")
#RNAseq.logFC <- merge(RNAseq.logFC, cuff2name, by.x = 0, by.y = "cuff_ids")

RNAseq.logFC.ruved <- Mm.RUVg.model[[1]][, c("logFC.N3vs0", "logFC.N5vs0",
                                             "logFC.N7vs0", "logFC.B5vs0",
                                             "logFC.R5vs0" )]

RNAseq.logFC.ruved <- merge(RNAseq.logFC.ruved, annot.frame,
                            by.x = 0, by.y = "ensembl_id")
#RNAseq.logFC.ruved <- merge(RNAseq.logFC.ruved, cuff2name,
#                            by.x = 0, by.y = "cuff_ids")

RNAseq.Array.logFC <- merge(RNAseq.logFC, Array.logFC,
                            by.x = "gene_name", by.y = 0) # NOTE!!!! gene_name or gene_names

RNAseq.Array.logFC <- merge(RNAseq.Array.logFC, RNAseq.logFC.ruved,
                            by = c("Row.names","gene_name"))

names(RNAseq.Array.logFC) <- gsub(".x", ".plain", names(RNAseq.Array.logFC))

names(RNAseq.Array.logFC) <- gsub(".y", ".ruved", names(RNAseq.Array.logFC))
## Very interesting...
cor(RNAseq.Array.logFC[,3:14])
cor(RNAseq.Array.logFC[,3:14], method="spearman")

pdf("figuresANDmanuscript/Array_vs_RNAseq_pairs.pdf")
#ggpairs(RNAseq.Array.logFC[, 3:14], alpha=0.2) + theme_bw()
ggpairs(RNAseq.Array.logFC[, 3:14], mapping=aes(alpha=0.2)) + theme_bw()
dev.off()

##Warning messages when script is run:
##2: In ggpairs(RNAseq.Array.logFC[, 3:14], alpha = 0.2) :
## Extra arguments: 'alpha' are being ignored.  
## If these are meant to be aesthetics, submit them using the 'mapping' variable within 
## ggpairs with ggplot2::aes or ggplot2::aes_string.
##
pdf("figuresANDmanuscript/Array144vsRNAseqN7.pdf")
ggplot(RNAseq.Array.logFC, aes(X144h, logFC.N7vs0.plain)) +
    geom_point(alpha=0.8) +
        stat_density2d(aes(alpha=..level.., fill=..level..), size=2,                                                                          geom="polygon") +
            scale_fill_gradient(low = "yellow", high = "red") +
                scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
                    stat_smooth() +
                        theme_bw()
dev.off()

pdf("figuresANDmanuscript/Array144vsRNAseqN7_ruved.pdf")
ggplot(RNAseq.Array.logFC, aes(X144h, logFC.N7vs0.ruved)) +
    geom_point(alpha=0.8) +
        stat_density2d(aes(alpha=..level.., fill=..level..), size=2,                                                                          geom="polygon") +
            scale_fill_gradient(low = "yellow", high = "red") +
                scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
                    stat_smooth() +
                        theme_bw()
dev.off()


