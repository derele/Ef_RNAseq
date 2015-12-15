## read the array data as done for Schmidt et al.  "Eimeria
## falciformis infection of the mouse caecum identifies opposing roles
## of IFNg-regulated host pathways for the parasite development"
## DOI:10.1038/mi.2013.115

library(mgug4122a.db)
library(Agi4x44PreProcess)
library(Biobase)
library(RSvgDevice)
library(GGally)

targets <- readTargets("/data/Mouse_arrays/2012_arrays/targets.txt", sep=";")

RG <- read.maimages(targets, path= "/data/Mouse_arrays/2012_arrays/", source = "agilent", 
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

biolrep <- targets$SampleNumber

design <- modelMatrix(targets, ref="Uninf")
cont.matrix <- makeContrasts("X24h" , "X144h", "X24h-X144h",
                             levels=design)

corfit.g <- duplicateCorrelation(MA.g, design=design, ndups=1, block=biolrep)

fit.g <- lmFit(MA.g, design=design, block=biolrep, cor=corfit.g$consensus)

fit2.g <-  contrasts.fit(fit.g, cont.matrix)
fit2.g <- eBayes(fit2.g)

Array.logFC <- topTable(fit2.g, n=400000)[, c("X24h", "X144h")]

RNAseq.logFC <- ALL.top.ALT[, c("logFC.N3vs0", "logFC.N5vs0",
                                "logFC.N7vs0", "logFC.B5vs0", "logFC.R5vs0" )]
RNAseq.logFC <- merge(RNAseq.logFC, id2name, by.x = 0, by.y = "gene_ids")

RNAseq.Array.logFC <- merge(RNAseq.logFC, Array.logFC, by.x = "gene_names", by.y = 0)


## Very interesting...
cor(RNAseq.Array.logFC[,3:9])

ggpairs(RNAseq.Array.logFC[, 3:9], alpha=0.2) + theme_bw()

ggplot(RNAseq.Array.logFC, aes(logFC.N7vs0, X144h)) +
    geom_point(alpha=0.8) +
        stat_density2d(aes(alpha=..level.., fill=..level..), size=2,                                                                          bins=10, geom="polygon") +
            scale_fill_gradient(low = "yellow", high = "red") +
                scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
                    stat_smooth() +
                        theme_bw()

## ... for a distinct subset of highly DE genes in RNAseq the DE was
## not found on Array
