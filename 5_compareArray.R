## RNAseq - Microarray Comparison
## read the array data as done for Schmidt et al.  "Eimeria
## falciformis infection of the mouse caecum identifies opposing roles
## of IFNg-regulated host pathways for the parasite development"
## DOI:10.1038/mi.2013.115

library(mgug4122a.db)
library(Agi4x44PreProcess) # !!!! NOT available for R 3.3.0 - NOT on Bioconductor anymore
library(Biobase)
library(RSvgDevice)
library(GGally)
library(ggplot2)

## read the fold-change data for mouse
Mm.DE.test <- read.table("output_data/Mm_DEtest.csv", sep=",")

## make it wider format to merge with array results (bleow)
Mm.DE.FC <- Mm.DE.test[, c("contrast", "gene", "logFC")]
Mm.DE.FC <- reshape(Mm.DE.FC, idvar = "gene", timevar = "contrast", direction="wide")

### getting all the annotation names to be able to merge based on them
annot.frame <- read.table("output_data/annotation_data.csv", sep=",")

targets <- readTargets("data/targets.txt", sep=";")

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

Array.logFC <- merge(annot.frame[,1:3], Array.logFC, 
                     by.x = "symbol", by.y = 0)

RNAseq.Array.logFC <- merge(Array.logFC[,c("ensembl_id", "X24h", "X144h")],
                            Mm.DE.FC, 
                            by.x = "ensembl_id", by.y = "gene")


RNAseq.Array.logFC <- unique(RNAseq.Array.logFC)

cor(RNAseq.Array.logFC[, -1], method="spearman")[,1:2]

cor.test(RNAseq.Array.logFC[, "logFC.N7vs0" ],
         RNAseq.Array.logFC[, "X144h",], method="spearman", exact = FALSE)

#pdf("figures/SI_1_Array144vsRNAseqN7.pdf", onefile = FALSE)
array.comarison <- ggplot(RNAseq.Array.logFC, aes(X144h, logFC.N7vs0)) +
  geom_point(alpha=0.8) +
  stat_density2d(aes(fill=..level..), size=2, geom="polygon") +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_alpha(range = c(0.00, 0.95), guide = FALSE) +
  stat_smooth() +
  theme_bw(20, size = 20) +
  xlab("Microarray data 6 dpi") +
  ylab("RNA-seq data 7dpi")
ggsave("Supplement/SI_5_QC_arrayCorrelation_MDS/SI_Array144vsRNAseqN7.svg", height = 12, width = 12, plot = array.comarison)
#dev.off()


