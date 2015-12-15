### ballgown for statistical DEG analysis of RNAseq data
## mouse and Eimeria falciformis as input


library(ballgown)
setwd("~/Ef_RNAseq")

## IMPORT DATA USING BALLGOWN IMPORT
if(!exists("rnaseq.bg")) {
    rnaseq.bg = ballgown(dataDir="tablemaker_060815", samplePattern = "*")
}

source("srcEH/functionsEfRNAseq.R")

########################

## basically design matrix, providing ballgown with the experimental
# gene_expression = gexpr(rnaseq.bg) # can remove???

add.pdata <- function (ballgown.obj){
    samples <- sub("^tablemakertophat_(.*)(_hiseq|_forw.fastq_paired)$",
                   "\\1", sampleNames(ballgown.obj))
    grouped <- sub("_rep\\d+$" , "", samples)
    infection <- ifelse(grepl("1stInf", grouped), "1st",
                        ifelse(grepl("2ndInf", grouped), "2nd", "NULLth"))
    rep <- sub(".*_(rep\\d+)$" , "\\1", samples) # (name here) is taken "\\1" here
    mouse.strain <- sub("^(.*?)_.*", "\\1", grouped)
    timepoint <- as.numeric(as.character(sub("^.*_(.*)dpi", "\\1", grouped)))
    ## variables to summarize other variables sith mutliple levels into two levels ##
    immune.status <- ifelse(mouse.strain %in% "Rag", "Rag", "competent")
    late.early <- ifelse(timepoint < 7 , "early", "late")
    ##                                       ##
    data.frame(samples, grouped, infection, rep, mouse.strain, timepoint, immune.status, late.early)
}

pData(rnaseq.bg) <- add.pdata(rnaseq.bg)

########################################
## Create mouse only object and remove sporozoite and oocyst samples from mouse data
mouse.bg <- subset(rnaseq.bg, grepl('^XLOC.*', geneIDs(rnaseq.bg)),
                   genomesubset=TRUE)
mouse.bg <- subset(mouse.bg, !is.na(pData(mouse.bg)$timepoint),
                   genomesubset=FALSE) # 

## Subset only infected samples 
inf.mouse.bg <- subset(mouse.bg, !pData(mouse.bg)$timepoint==0,
                       genomesubset=FALSE)

## Subset uninfected samples 
uinf.mouse.bg <- subset(mouse.bg, pData(mouse.bg)$timepoint==0,
                        genomesubset=FALSE)

## Subset only the timecourse experiment: Infected in NMRI
NMRI.I.mouse.bg <- subset(inf.mouse.bg,
                          pData(inf.mouse.bg)$mouse.strain%in%"NMRI",
                          genomesubset=FALSE)

## Subset for dpi5 for which we have the other mouse strains 
dpi5.mouse.bg <- subset(mouse.bg,
                        pData(mouse.bg)$timepoint==5,
                        genomesubset=FALSE)

## Create Eimeria only object and remove day zero
Ef.bg <- subset(rnaseq.bg,
                grepl('^EfaB.*', geneIDs(rnaseq.bg)),
                genomesubset=TRUE)
Ef.bg <- subset(Ef.bg,
                pData(Ef.bg)$timepoint != 0 | is.na(pData(Ef.bg)$timepoint),
                genomesubset=FALSE)


##################################
##  M     O      U      S     E ##
##################################

###############################            TIME              ##########################

## time differences in all
stat_results_time = stattest(mouse.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_time = stat_results_time[order(stat_results_time$qval),]
##  significant stuff...  

## Only looking at infected mouse data
stat_results_timeI <- stattest(inf.mouse.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_timeI <- stat_results_timeI[order(stat_results_timeI$qval),]
## all significance gone down the drains ... as the significance is
## for infected / non infected

## Looking at infected mouse data but without the potentiall unnecessary strain covariate
stat_results_timeNI <- stattest(NMRI.I.mouse.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection"), 
    timecourse = F) 
stat_results_timeNI <- stat_results_timeNI[order(stat_results_timeNI$qval),]
## pheewww, further down the drains ... as the significance is for
## infected / non infected and the non dpi5 data gave a) higher rep
## sizes OR b) just messed things up and gave us unjustified confidence

## Only looking at infected mouse data but summarizing time to early late comparison
stat_results_timeLE <- stattest(inf.mouse.bg, feature='transcript',
    meas='cov', 
    covariate = 'late.early', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_timeLE <- stat_results_timeLE[order(stat_results_timeLE$qval),]
##  still no significance 

## NO SIGNIFICANT DIFFERENCES in mice EARLY/LATE IN infections.


###########################                  STRAINS             #################
## contrasting immune vs. non-immune in all
## compare immune.status with mouse.strain output later...

stat_results_strain = stattest(mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection", "timepoint"), 
                               timecourse = F) 
stat_results_strain = stat_results_strain[order(stat_results_strain$qval),]
## large set with high significance

## contrasting immune vs. non-immune in only day 5 
stat_results_strain5 = stattest(dpi5.mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection"), 
                               timecourse = F) 
stat_results_strain5 = stat_results_strain5[order(stat_results_strain5$qval),]
## smaller set with significance but this can be as we use power...

## running on only infected samples 
stat_results_strainI = stattest(inf.mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection", "timepoint"), 
                               timecourse = F) 
stat_results_strainI = stat_results_strainI[order(stat_results_strainI$qval),]
## smaller set with significance but this can be as we use power...

## running on only UNinfected samples 
stat_results_strainUI = stattest(uinf.mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection"), 
                               timecourse = F) 
stat_results_strainUI = stat_results_strainUI[order(stat_results_strainUI$qval),]
## agin smaller set with significance but this can be as we use FURTHER power...

## Lets cross-check these: We are more interesed in genes showing
## difference in the strains when infected rather than a plain
## difference when uninfected 

## How many of the top 50/100 in day5 are in top X uninfected differces
table (head(stat_results_strain5$id, 50) %in% head(stat_results_strainUI$id, 500))
table (head(stat_results_strain5$id, 50) %in% head(stat_results_strainUI$id, 1000))
table (head(stat_results_strain5$id, 100) %in% head(stat_results_strainUI$id, 500))
table (head(stat_results_strain5$id, 100) %in% head(stat_results_strainUI$id, 1000))
## Okay
## What about the different genes in infected samples
table (head(stat_results_strainI$id, 50) %in% head(stat_results_strainUI$id, 500))
table (head(stat_results_strainI$id, 50) %in% head(stat_results_strainUI$id, 1000))
table (head(stat_results_strainI$id, 100) %in% head(stat_results_strainUI$id, 500))
table (head(stat_results_strainI$id, 100) %in% head(stat_results_strainUI$id, 1000))
## So this is Okay. We see differences in mouse immunocompetent vs
## incompetent mouse strains upon infection


################## first vs. second infection ### #############################
stat_results_chal = stattest(mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'infection', 
                               adjustvars = c("timepoint", "immune.status"), 
                               timecourse = F) 
stat_results_chal = stat_results_chal[order(stat_results_chal$qval),]


## only in infected mice
stat_results_chalI = stattest(inf.mouse.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'infection', 
                               adjustvars = c("timepoint", "immune.status"), 
                               timecourse = F)
stat_results_chalI = stat_results_chalI[order(stat_results_chalI$qval),]
## Hey there is someting ant its not a result from infected / uninfected

## have these transcripts something in common with the RAG / WT differences
## -> adaptive immune system ...
## How many of the top 50/100 in day5 are in top X uninfected differces
intersect (head(stat_results_strain5$id, 50) , head(stat_results_chalI$id, 50))
intersect (head(stat_results_strain5$id, 100) , head(stat_results_chalI$id, 100))
intersect (head(stat_results_strain5$id, 500) , head(stat_results_chalI$id, 500))
intersect (head(stat_results_strain5$id, 1000) , head(stat_results_chalI$id, 1000))
## ... not really ... or do they: Very basic stats ;-).

## if yes maybe...
immuneGenes <- intersect (head(stat_results_strain5$id, 1000) , head(stat_results_chalI$id, 1000))
immuneGenes <- as.character(mouse.bg@indexes$t2g[mouse.bg@indexes$t2g$t_id%in%immuneGenes,
                                  "g_id"])

###########################################

## Get DE genes
## chose q-value threshold, multiple correction adjusted p-val
## check: table(stat_results_inf$qval<0.05)

select.from.stats.results <- function (stats.results, bg.obj, qval=0.05){
    diff_id <- stats.results[which(stats.results$qval<qval), "id"]
    diff_id <- as.character(diff_id)
    diff_genes <- bg.obj@indexes$t2g[bg.obj@indexes$t2g$t_id%in%diff_id,
                                     "g_id"]
    return(as.character(diff_genes))
}

## time_genes <- select.from.stats.results(stat_results_time, mouse.bg)
## well ther are none really...

strain_genes5 <- select.from.stats.results(stat_results_strain5, mouse.bg, qval = 0.05)
strain_genesA <- select.from.stats.results(stat_results_strain, mouse.bg, qval = 0.05)
strain_genesI <- select.from.stats.results(stat_results_strainI, mouse.bg, qval = 0.05)

### hey try 0.01 and compare to 0.05 what a wired distribution of
### qvalues!

chal_genes <- select.from.stats.results(stat_results_chal, mouse.bg)
chal_genesI <- select.from.stats.results(stat_results_chalI, mouse.bg)

### hey try 0.01 and compare to 0.05 what even more wired distribution
### of qvalues!




