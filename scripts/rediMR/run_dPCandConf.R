#dietPC adjustments

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


# load required packages
suppressMessages(silent <- lapply(c("tidyverse", "data.table", "parallel", 
                                    "paletteer", "RColorBrewer", "ggpubr"),  
                                  library, character.only = T))

# Run pantry.R file
source("../scripts/pantry.R")


#############################################
## ~~~~~~ Set up & input parameters ~~~~~~ ##
#############################################

# ==========================
## Load command arguments
# ==========================

# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1] #diet phenotype
tag <- args[2] #version "tag"
covarGroup <- args[3]


pheno_tag <- paste0(pheno, "_", tag)
pheno_tag_covarGroup <- paste0(pheno_tag, "_", group)

redimrDir <- "../data/processed/rediMR"
ssInput <- paste0(redimrDir, "/", pheno_tag, "_ssInput.csv")
datInput <- paste0(redimrDir, "/", pheno_tag, "_datInput.rda") 


## Inputs & parameters
inputs <- paste0("\n compare_rediMR Inputs: \n -pheno ", pheno, "\n -tag ", tag,  
                 "\n -group ", covarGroup)
files <- paste0("\n -ssInput: ", ssInput, "\n -datInput: ", datInput, 
                 "\n -redimrDir ", redimrDir)

cat(inputs)
cat(files)


# ==========================================
## Load data input files
# ==========================================

## Input files
ss<-fread(ssInput) %>% filter(LOCI == 1)
dat<-readRDS(datInput)


## covariate sets
covarSetGroups[[covarGroup]]$Sets

cat(paste0("Comparing the following covariate groups:"), 
    paste0("\n - ", covarSetGroups[[covarGroup]]$Sets), "\n")


####################################################################################################
## Examine impact of adjusting for additional lifestyle/SES covariates on SNP-diet associations
####################################################################################################

snps <- names(dat %>% select(starts_with("rs") | starts_with("snp")))

# ===================================================================
## Adjusting SNP-food associations for ALL dietPCs AND confounders 
# ===================================================================

# Dot plot of pctBchange when adjusting for ALL dietPCs PLUS confounders =============

BchangeEachConf <- lapply(24:30, function(set){
  do.call(rbind.data.frame, lapply(snps, function(snp) {
    pctBchange.fun(pheno, snp, adjCovar = covarSets[[ covarSetGroups$dietEachConf$Sets[[set]] ]]$Covars, 
                   covarName = covarSets[[ covarSetGroups$dietEachConf$Sets[[set]] ]]$Label, data=dat)
  }) )
})

BchangeEachConf %>% saveRDS(paste0("../data/processed/rediMR/", pheno, "_BchangeEachConf.rda"))


# Dot plot of pctBchange when adjusting for ALL dietPCs PLUS ADDITIONAL confounders ============

BchangeAddConf <- lapply(24:30, function(set){
  do.call(rbind.data.frame, mclapply(snps, function(snp) {
    pctBchange.fun(pheno, snp, adjCovar = covarSets[[ covarSetGroups$dietAddConf$Sets[[set]] ]]$Covars, 
                   covarName = covarSets[[ covarSetGroups$dietAddConf$Sets[[set]] ]]$Label, data=dat)
  }) )
})

BchangeAddConf %>% saveRDS(paste0("../data/processed/rediMR/", pheno, "_BchangeAddConf.rda"))



#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


## EOF



