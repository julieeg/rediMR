# ReDiMR - Step 1: SNP Refinement using Covariate Adjustment

# load required packages
lapply(c("tidyverse", "data.table", "paletteer", "RColorBrewer", "ggpubr"),  
       library, character.only = TRUE)


# load pantry file with stored parameters, basic functions & ggplot templates
source("../scripts/pantry.R")


#############################################
## ~~~~~~ Set up & input parameters ~~~~~~ ##
#############################################

# ==========================
## Load command arguments
# ==========================

args <- commandArgs(trailingOnly = TRUE)
exposure = args[1] # "oilyfish_QT"
outcome = args[2] # "tg" 
sumstats = args[3] # ../data/sumstats/standardized_shared_from_KSJ/std_mrdat_oilyfish_GCST90239664_TG_Graham_GRCh37.csv
covarset = args[4] #"confounders1"
tag = args[5] # "vCole" 
saveDir= args[6] #../data/processed/rediMR/vCole/oilyfish_QT_tg

phenofile = "../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits_09292024.txt"
genofile = paste0(saveDir, "/", exposure, "_", outcome, "_snpsInput.raw")
covars_gwas="age sex gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10" #args[3]

savePref = paste0(exposure, "_", outcome, "_", covarset)


## Inputs & parameters
inputs <- paste0("\n ReDiMR Step 1 Inputs: \n -exposure ", exposure, 
    "\n -outcome ", outcome,
    "\n -sumstats ", sumstats,
    "\n -covars_gwas ", covars_gwas,
    "\n -covarset ", covarset, 
    "\n -phenofile ", phenofile, 
    "\n -genofile ", genofile, 
    "\n -saveDir ", saveDir,
    "\n -savePref ", savePref,
    "\n -tag ", tag,
    "\n"
    )

cat(paste0(inputs))


# ==========================================
## Load data input files
# ==========================================

# genotype file with dosage
dose <- fread(genofile)
dose_id <- dose %>% rename(id=FID) %>% select(-c(IID, MAT, PAT, SEX, PHENOTYPE )) %>%
  rename_with(., ~gsub(":",".", .)) %>%
  rename_with(., ~ifelse(startsWith(., "rs") | . == "id", ., paste0("id",.)))

# Rename SNPs with chr:pos_A1_A2 naming structure
snps <- names(dose_id %>% select(-id))

# phenotype file with covariates 
phenos_id <- fread(phenofile)

# Define covariate set
covarSet <- covarSets[[covarset]]

# ==================================
## Build dietPCs without exposure 
# ==================================
 
if(any(startsWith(covarSet$Covars, "dietPC"))) {
  
  if(exposure == "oilyfish_QT") { 
    exclude <- "oily_fish"
  } else if(exposure == "bread_type_BIN") { 
    exclude <- c("bread_type_white_vs_brown_or_whole", "bread_intake")
  } else if (exposure == "alch_glasspermonth_QT") {
    exclude <- "none"
  }
  
  dietpcs = derive_dietPCs(exposure, exclude, data=phenos_id)
  
  # Extract dietPC scores
  dietPCs.scores <- dietpcs$scores
  dietPCs.loadings <- dietpcs$pcs

  # Add to phenotype file
  phenos_id <- phenos_id %>% left_join(dietPCs.scores, by="id")

  dietPCs.loadings %>% saveRDS(paste0(saveDir, "/", savePref, "_dietPCloadings_", tag, ".rda"))
  
}


# ==================================
## Compile dataset
# ==================================

phenos_id %>% 
  select(id, 
         all_of(exposure),
         all_of(strsplit(covars_gwas, split=" ")[[1]]),
         all_of(covarSet$Covars)
)

dat <- left_join(phenos_id, dose_id, by="id")

cat("\nCovariates for gwas:", covars_gwas, 
    "\nCovariats for adjustment:", paste0("\n   - ",
                                          covarSet$Names, " = ", covarSet$Covars)
)

##########################################################
## STEP 1: SNP Refinement based on covariate adjutment  ## 
##########################################################

cat("\n Starting Step 1: SNP Refinement ... \n ")

# ================================================
## Adjust for ALL covariates --> SNP Refinement  
# ================================================

cat("\n Calculating %change in Beta when adjusting for ALL covariates, relative to a BASE model")

# For eacn SNP, tabulate pctBchange when adjusting for ALL covariates
covars_gwas.vars <- strsplit(covars_gwas, " ")[[1]]
covars_adjust.vars <- covarSet$Covars #strsplit(covars_adjust, " ")[[1]]

tab_bchange_all <- do.call(rbind.data.frame, lapply(snps, function(snp) {
  pctBchange.fun(pheno=exposure, snp, 
                 adjCovar=covars_adjust.vars, baseCovars=covars_gwas.vars,
                 covarName = "All_Covariates", data=dat) } ))

# Write results to csv
head(tab_bchange_all)


# =======================================================
## Adjusting for EACH covariate --> Covariate Influence  
# =======================================================

cat("Calculating pctBchange when adjusting for EACH covariate ... \n")

# For each snp, tabulate pctBchange when adjusting for EACH covariate
tab_bchange_each <- do.call(rbind.data.frame, lapply(snps, function(snp) { 
  do.call(rbind.data.frame, lapply(covars_adjust.vars, function(covar) {
    pctBchange.fun(pheno=exposure, snp, 
                   adjCovar=covar, 
                   covarName=covar, data=dat) } ))  } ))

# Write results to csv
head(tab_bchange_each)



# =========================
## Compile & save as .csv
# =========================

## Compile rediMRdat for downstream analyses
rbind.data.frame(tab_bchange_all, tab_bchange_each) %>%
  write.csv(paste0(saveDir, "/", savePref, "_bchange_full_",tag,".csv"))


##########################
## STEP 1 COMPLETED !  ##
##########################

cat(paste0(
"You completed ReDiMR Step 1: SNP Refinement for ", exposure, " & ", outcome, "summary statistics, 
with covariate set ", covarset, ", \n", "Now proceeding to Step 2: Two-Sample MR ... "))


#EOF

