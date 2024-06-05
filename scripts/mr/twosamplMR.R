# ReDiMR
# Last updated: June 3, 2023


#################################
## Set up & assign parameters  ##
#################################

# =======================
##  Set up
# =======================

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
tag <- args[2]
#outcome_id <- args[3]

pheno_tag <- paste0(pheno, "_", tag)

rediMR_inputPref=paste0("../data/processed/rediMR/", pheno_tag)
ssInput <- paste0(rediMR_inputPref, "_ssInput.csv") #args[3]
datInput <- paste0(rediMR_inputPref, "_datInput.rda")  #args[4] 
outDir <- "../data/processed/mr" #args[6]


# load basic functions
source("../scripts/basic_functions.R")


#################################################################
##  Load input files ; Write functions to calculate %B change  ##
#################################################################

# ==========================================
## Load data input files
# ==========================================

# Load summary stats data
ss<-fread(ssInput) %>% filter(LOCI == 1)

# Load phenotype/dosage data
dat<-readRDS(datInput)

# B pct change from all covariate adjustment
tabBchangeAllCov <- fread(paste0(rediMR_inputPref, "_tabBchangeAllCov.csv"))



##########################################
##  STEP 2: Set up & assign parameters  ##
##########################################

cat("Starting Step 2: Two-Sample Mendelian Randomization ... \n ")


remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github('MRCIEU/ieugwasr')
library(TwoSampleMR) ; library(ieugwasr)


# establish token for ieugwas access
opengwas_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqdWxpZS5nZXJ2aXNAdHVmdHMuZWR1IiwiaWF0IjoxNzE3MTc5MzkxLCJleHAiOjE3MTgzODg5OTF9.Uu_WDLV7LmYN9ZYoEnVr_pOEUu9ggPHLH3gdUNJ4CcyBSG4VQ2eB-prv8SzaCINn8-pEExBBSnrUx6lHS5324g9qvQWV52cVZuwBEXjlvqxKz99Lq6AIRsUWaTxAZ38JWIiGF6M2bElty9G0NE4WPIWi76OTE9Ad4iv5SQZvFWO-4qwxzxudcWnsNasvd6Y12HSaTnMjDenR2DU8guF_O_bFL3yRxOjCKZT6WKfplmRZXJM8HADamK7h-FUUp7RMavgvLeLWjUhVT1pLVL3BvHSvOiH49lRl9FBG_KAuO6nnewWZCVh2eAof9L2Ud2kec4o_PIF8VsoAl8bwe-aUvA"


# ======================================================
## Build exposure_dat file for TwoSampleMR
# ======================================================

# -tabBchangeAllCov for refined/included SNP subsets
# -loci.afreq for effect allele frequency
# -ssInput.csv for [unrefined] gwas summary stats
afreq <- fread(paste0("../data/processed/rediMR/", pheno_tag, "_loci.afreq")) %>%
  #mutate(SNP=ifelse(!startsWith(ID, "rs"), gsub(":", ".", paste0("snp", ID)), ID)) %>%
  mutate(REF_FREQS = 1-ALT_FREQS) %>% 
  select(SNP=ID, ALT_FREQS, REF_FREQS, REF_FREQ_ALLELE=REF, ALT_FREQ_ALLELE=ALT) 


exposure_dat <- ss %>% rowwise() %>% 
  mutate(refined.set = tabBchangeAllCov$RefinedSet[which(startsWith(tabBchangeAllCov$snp, SNP))]) %>%
  ungroup() %>% 
  mutate(SNP=gsub("snp", "", gsub("[.]", ":", SNP))) %>%
  left_join(afreq, by = "SNP") %>% 
  mutate(exposure = paste0(pheno), id.exposure=paste0("ukb-", pheno),
         effect_allele.exposure = ifelse(BETA>0, EA, NEA),
         other_allele.exposure = ifelse(BETA>0, NEA, EA),
         beta.exposure=abs(BETA), se.exposure=SE) %>%
  mutate(eaf.exposure=ifelse(effect_allele.exposure==ALT_FREQ_ALLELE, ALT_FREQS, REF_FREQS)) %>%
  select(SNP, exposure, id.exposure, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, refined.set)


# compile as list by snp refinement set
exposure_dat.l <- lapply(list(c(1,0),1,0), function(set) exposure_dat %>% filter(refined.set %in% set))
snp_sets <- c("All", "Refined", "Unrefined")
exposure_dat.l <- lapply(1:3, function(set) exposure_dat.l[[set]] %>% mutate(snp_set=snp_sets[set]))



# ======================================
## Compile outcomes from ieuGWAS
# ======================================

outcomes <- list(
  cbind.data.frame(outcome_var="colorectal_cancer", id.outcome="finn-b-C3_COLORECTAL", type="or"),
  cbind.data.frame(outcome_var="ischemic_stroke", id.outcome="finn-b-I9_STR_EXH_EXNONE", type="or"),
  cbind.data.frame(outcome_var="depression", id.outcome="finn-b-F5_DEPRESSIO", type="or"),
  cbind.data.frame(outcome_var="t2d", id.outcome="ebi-a-GCST006867", type="or"),
  cbind.data.frame(outcome_var="crp", id.outcome="prot-a-670", type="beta"),
  cbind.data.frame(outcome_var="tg", id.outcome="ebi-a-GCST002216", type="beta"),
  cbind.data.frame(outcome_var="dbp", id.outcome="ebi-a-GCST90018952", type="beta"),
  cbind.data.frame(outcome_var="glucose", id.outcome="ieu-b-113", type="beta")
)

outcome_dat.l <- lapply(outcomes, function(out) {
  extract_outcome_data(snps=exposure_dat$SNP, outcome = out$id.outcome, rsq=0.5, opengwas_jwt = opengwas_jwt) %>%
    mutate(type=rep(out$type, nrow(.)))
}) ; names(outcome_dat.l) <- outcomes

nOutcomes=length(outcome_dat.l)

# Add indictor to outcome.l for whether estimates should be betas/ors
names(outcome_dat.l) <- c(sapply(outcomes, function(x) x$outcome_var))



##################
##    Run MR    ##
##################

# Write function to run MR for an exposure_dat & outcome_dat
## Run over list of SNPs ------------------------------------------------

mr_results.l <- list()
run_mr <- function(exposure_dat, outcome_dat) {
  # run MR
  mr_dat=harmonise_data(exposure_dat, outcome_dat) %>%
    mutate(Fstat=(beta.exposure^2)/((se.exposure)^2),.after=se.exposure)
  mr_res=mr(mr_dat) 
  
  summary=mr_res %>% 
    select(id.exposure, id.outcome, outcome, method, nsnp, b, se, pval) %>% 
    mutate(snp_set=exposure_dat$snp_set[1], .before="id.exposure") %>%
    mutate(id.exposure=paste0(id.exposure, "-", snp_set, "_snps")) %>%
    mutate(outcome=gsub("[||].*","",outcome), .after=id.outcome) %>%
    mutate(lci95=-1.96*se, uci95=b+1.96*se, .after="se")
  
  # If nsnp>1, calculate heterogeneity & pleiotropy
  if(summary$nsnp[1] >1) {
    het=mr_heterogeneity(mr_dat) ; pleio=mr_pleiotropy_test(mr_dat)
    summary <- summary %>% left_join(het %>% select(method, Q, Q_df, Q_pval)) %>% 
      left_join(pleio %>% mutate(method = "MR Egger") %>% 
                  select(method, egger_intercept, egger_se=se, egger_pval=pval), by = "method") %>%
      mutate(method=factor(method, levels=c("MR Egger", "Inverse variance weighted", "Weighted median", "Weighted mode", "Simple mode"))) %>%
      arrange(method)
    } else if (summary$nsnp[1]<=1) { summary <- summary %>% bind_cols(
      as.data.frame(matrix(NA, 1, 6, dim=list(NULL, c("Q", "Q_df", "Q_pval", "egger_intercept", "egger_se", "egger_pval")))))
    }
  
  #Transform beta to OR if specified; Rename b/or as 'estimate' for consistency
  if(outcome_dat$type[1] !="beta") {
    summary <- summary %>% select(-c(b,se,lci95,uci95,pval)) %>% mutate(
      generate_odds_ratios(mr_res) %>% select(estimate=or, se, lci95=or_lci95, uci95=or_uci95, pval), .after=nsnp) %>% 
      mutate(type=rep(outcome_dat$type[1], nrow(.)), .after=outcome)
  } else { summary <- summary %>% rename("estimate"=b) %>% mutate(type="beta", .after=outcome) }
  
  # Compile MR results as list
  mr_results <- list(
    mr_dat=mr_dat, mr_res=mr_res, mr_summary=summary
  )
  return(mr_results)
}


# Run MR over list of outcomes for each exposure snp set
mr_results.l <- lapply(outcome_dat.l, function(outcome) {
  sets.l <- lapply(exposure_dat.l, function(exposure) {
    run_mr(exposure, outcome) }) ; names(sets.l)<-c("All", "Refined", "Unrefined")
    return(sets.l)})

lapply(as.list(1:nOutcomes), function(i) {
  saveRDS(mr_results.l[[i]], file = paste0(outDir, "/", pheno_tag, "_MRsummary_", names(outcome_dat.l)[i], ".rda"))
})



###########################################################
##  Run MR Sensitivity Analysis - Validity vs. Strength  ##
###########################################################






###############################################################################
                            ## ~~~ Step 2 Completed !! ~~~~ ##
cat(paste0("\n Congratulatulations!!\n",
           " You completed ReDiMR Step 2: Two-Sample MR for ", pheno, "!"))
###############################################################################


#EOF


