# ReDiMR
# Last updated: June 12, 2024


#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#


# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# load pantry file with stored parameters, basic functions & ggplot templates
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

pheno_tag <- paste0(pheno, "_", tag)

covar <- args[3]
ssInput <- args[4] #paste0("../data/processed/rediMR/", pheno_tag, "_ssInput.csv")
datInput <- args[5] #paste0("../data/processed/rediMR/", pheno_tag, "_datInput.rda") 
pctBthold <- args[6] #20
redimrDir <- args[7]
gwasCovars <- args[8]

pheno_tag_covarSet <- paste0(pheno_tag, "_", covar)


## Inputs & parameters
inputs <- paste0("\n ReDiMR Inputs: \n -pheno ", pheno, "\n -tag ", tag,  
    "\n -ssInput ", ssInput, "\n -datInput ", datInput, 
    "\n -gwasCovars ", gwasCovars,
    "\n -covarSet ", covar, "\n -pctBthold ", pctBthold, 
    "\n -redimrDir ", redimrDir)

cat(paste0(inputs))


# ==========================================
## Load data input files
# ==========================================

# Load summary stats data
ss<-fread(ssInput) %>% filter(LOCI == 1)

# Load phenotype/dosage data
dat<-readRDS(datInput)

# list of snps
snps <- names(dat %>% select(contains(ss$SNP)))

# covariate sets
covarSet <- covarSets[[covar]]

if(length(covarSet$Names)>1) {
  cat("\n Covariates for adjustment:", covarSet$Label, "Set \n", 
      paste("-", covarSet$Names, "\n") )
}
cat("\n Covariates for adjustment:", covarSet$Label)




#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#


##########################################################
## STEP 1: SNP Refinement based on covariate adjutment  ## 
##########################################################

cat("\n Starting Step 1: SNP Refinement ... \n ")


# ================================================
## Adjust for ALL covariates --> SNP Refinement  
# ================================================

cat("\n Calculating pctBchange when adjusting for ALL covariates, relative to a BASE model
    with:", gwasCovarsBase, "... \n ")


# For eacn SNP, tabulate pctBchange when adjusting for ALL covariates
allCovars <- covarSet$Covars

tabBchangeAllCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) {
  pctBchange.fun(pheno=pheno, snp=snp, 
                 adjCovar=covarSet$Covars,
                 baseCovars=gwasCovars,
                 covarName = covarSet$Label, data=dat) }, mc.cores = 8 ))
head(tabBchangeAllCov)
tabBchangeAllCov <- tabBchangeAllCov %>% 
  mutate(RefinedSet = ifelse(abs(B_pctChange) < as.numeric(pctBthold),1,0))


# Write results to csv
write.csv(tabBchangeAllCov, file = paste0(redimrDir, "/", pheno_tag_covarSet, "_tabBchangeAllCov.csv"), row.names=T)


cat (paste0("DONE: Results written to ", paste0(redimrDir, "/", pheno_tag_covarSet, "_tabBchangeAllCov.csv")))
head(tabBchangeAllCov)



# =======================================================
## Adjusting for EACH covariate --> Covariate Influence  
# =======================================================

cat("Calculating pctBchange when adjusting for EACH covariate ... \n")

# For each snp, tabulate pctBchange when adjusting for EACH covariate
tabBchangeByCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) { 
  do.call(rbind.data.frame, mclapply(1:length(covarSet$Covars), function(i) {
    pctBchange.fun(pheno=pheno, snp=snp, 
                   adjCovar=covarSet$Covars[[i]], 
                   covarName=covarSet$Names[[i]], data=dat) }, mc.cores = 8))  }, mc.cores = 8))

# Write results to csv
fwrite(tabBchangeByCov, file = paste0(redimrDir, "/", pheno_tag_covarSet, "_tabBchangeByCov.csv"))


cat (paste0("DONE: Results written to ", paste0(redimrDir, "/", pheno_tag_covarSet, "_tabBchangeByCov.csv")))
head(tabBchangeByCov)



# =======================================================
## Make SNPset files for downstream analysis
# =======================================================

snpset <- tabBchangeAllCov %>% select(snp, RefinedSet) 

# rename SNPs to chr:pos format
snpset <- snpset %>% 
  mutate(ID=gsub("snp", "", gsub("[.]", ":", gsub("_[^_]*$", "", snp)))) %>%
  rename(SNP=snp)

snpset_ss <- ss %>% 
  mutate(ID=gsub("snp", "", gsub("[.]", ":", SNP))) %>%
  left_join(snpset, by = "ID")


# make snp sets
all <- snpset_ss$ID
refined <- (snpset_ss %>% filter(RefinedSet==1))$ID
unrefined <- (snpset_ss %>% filter(RefinedSet==0))$ID


## PRS input files ==========================

# All
snpset_ss %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag_covarSet, "_all_prsInput"))

# Refined
snpset_ss %>% 
  filter(RefinedSet == 1) %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag_covarSet, "_ref_prsInput"))

# Unrefined
snpset_ss %>% 
  filter(RefinedSet == 0) %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag_covarSet, "_unref_prsInput"))



##########################
## STEP 1 COMPLETED !!  ##
##########################

cat(paste0("
Congratulatulations!!\n",
"You completed ReDiMR Step 1: SNP Refinement for ", pheno, " with covariate set ",
covarSet$Label, " ! \n", "Now proceeding to Step 2: Two-Sample MR ... "))

################################################################################

#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#


#EOF

