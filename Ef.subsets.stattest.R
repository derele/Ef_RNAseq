## Subset EIMERIA transcripts for statistical tests, perform tests, and
## to select.from.stats.results function in annotations script


## Subset only infected samples not necessary:
## day 0 is excluded in Eimeria set

## Subset only the timecourse experiment: Infected in NMRI
NMRI.I.Ef.bg <- subset(Ef.bg,
                          pData(Ef.bg)$mouse.strain%in%"NMRI",
                          genomesubset=FALSE)

## Subset for dpi5 for which we have the other mouse strains 

########## problem #################

dpi5.Ef.bg <- subset(Ef.bg,
                        pData(Ef.bg)$timepoint==5,
                        genomesubset=FALSE)


################	 TIME		##############################
## time differences in all
stat_results_time = stattest(Ef.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_time = stat_results_time[order(stat_results_time$qval),]
##  significant stuff...  

## Only looking at infected samples
stat_results_timeI <- stattest(Ef.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_timeI <- stat_results_timeI[order(stat_results_timeI$qval),]
## Comment on significance in this set: 

## Looking at infected samples but without the potentially unnecessary strain
## covariate
stat_results_timeNI <- stattest(NMRI.I.Ef.bg, feature='transcript',
    meas='cov', 
    covariate = 'timepoint', 
    adjustvars = c("infection"), 
    timecourse = F) 
stat_results_timeNI <- stat_results_timeNI[order(stat_results_timeNI$qval),]
## Comment on significance in this set: 

## Only looking at infected mouse data but summarizing time to early late comparison
stat_results_timeLE <- stattest(Ef.bg, feature='transcript',
    meas='cov', 
    covariate = 'late.early', 
    adjustvars = c("infection", "immune.status"), 
    timecourse = F) 
stat_results_timeLE <- stat_results_timeLE[order(stat_results_timeLE$qval),]
## Comment on significance in this set: 

## General comment on time subsets for EIMERIA: 
##

###########################                  STRAINS             #################
## contrasting immune vs. non-immune in all
## compare immune.status with mouse.strain output later...

stat_results_strain = stattest(Ef.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection", "timepoint"), 
                               timecourse = F) 
stat_results_strain = stat_results_strain[order(stat_results_strain$qval),]
## Comment on significance in this set: 

## contrasting immune vs. non-immune in only day 5 
stat_results_strain5 = stattest(dpi5.Ef.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection"), 
                               timecourse = F) 
stat_results_strain5 = stat_results_strain5[order(stat_results_strain5$qval),]
## Comment on significance in this set: 

## running on only infected samples 
stat_results_strainI = stattest(Ef.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'immune.status', 
                               adjustvars = c("infection", "timepoint"), 
                               timecourse = F) 
stat_results_strainI = stat_results_strainI[order(stat_results_strainI$qval),]
## Comment on significance in this set: 

## How many of the top 50/100 in day5 are in top X uninfected differces
#table (head(stat_results_strain5$id, 50) %in% head(stat_results_strainUI$id, 500))
#table (head(stat_results_strain5$id, 50) %in% head(stat_results_strainUI$id, 1000))
#table (head(stat_results_strain5$id, 100) %in% head(stat_results_strainUI$id, 500))
#table (head(stat_results_strain5$id, 100) %in% head(stat_results_strainUI$id, 1000))
### Okay
## What about the different genes in infected samples
#table (head(stat_results_strainI$id, 50) %in% head(stat_results_strainUI$id, 500))
#table (head(stat_results_strainI$id, 50) %in% head(stat_results_strainUI$id, 1000))
#table (head(stat_results_strainI$id, 100) %in% head(stat_results_strainUI$id, 500))
#table (head(stat_results_strainI$id, 100) %in% head(stat_results_strainUI$id, 1000))
## So this is Okay. We see differences in mouse immunocompetent vs
## incompetent mouse strains upon infection


################## first vs. second infection ### #############################
## only in infected mice
stat_results_chalI = stattest(Ef.bg, feature='transcript',
                               meas='cov', 
                               covariate = 'infection', 
                               adjustvars = c("timepoint", "immune.status"), 
                               timecourse = F)
stat_results_chalI = stat_results_chalI[order(stat_results_chalI$qval),]
## Comment on significance in this set: 

## have these transcripts something in common with the RAG / WT differences
## -> adaptive immune system ...
## How many of the top 50/100 in day5 are in top X uninfected differces
#intersect (head(stat_results_strain5$id, 50) , head(stat_results_chalI$id, 50))
#intersect (head(stat_results_strain5$id, 100) , head(stat_results_chalI$id, 100))
#intersect (head(stat_results_strain5$id, 500) , head(stat_results_chalI$id, 500))
#intersect (head(stat_results_strain5$id, 1000) , head(stat_results_chalI$id, 1000))
## ... not really ... or do they: Very basic stats ;-).

## if yes maybe...
immuneGenes <- intersect (head(stat_results_strain5$id, 1000) , head(stat_results_chalI$id, 1000))
immuneGenes <- as.character(Ef.bg@indexes$t2g[Ef.bg@indexes$t2g$t_id%in%immuneGenes,
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

## time_genes <- select.from.stats.results(stat_results_time, Ef.bg)
## well ther are none really...

strain_genes5 <- select.from.stats.results(stat_results_strain5, Ef.bg, qval = 0.05)
strain_genesA <- select.from.stats.results(stat_results_strain, Ef.bg, qval = 0.05)
strain_genesI <- select.from.stats.results(stat_results_strainI, Ef.bg, qval = 0.05)

### hey try 0.01 and compare to 0.05 what a wired distribution of
### qvalues!

chal_genes <- select.from.stats.results(stat_results_chal, Ef.bg)
chal_genesI <- select.from.stats.results(stat_results_chalI, Ef.bg)

### hey try 0.01 and compare to 0.05 what even more wired distribution
### of qvalues!
