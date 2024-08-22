#MR excerpt from rediMR.R


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

covarSet <- args[3]
outcome <- args[4]

ssInput <- args[5] #paste0("../data/processed/rediMR/", pheno_tag, "_ssInput.csv")
datInput <- args[6] #paste0("../data/processed/rediMR/", pheno_tag, "_datInput.rda") 
pctBthold <- args[7] #20
mrDir <- args[8]

pheno_tag_covarSet <- paste0(pheno_tag, "_", covarSet)
tabBchangeAllCovInput <- paste0(redimrDir, "/", pheno_tag_covarSet, "_tabBchangeAllCov.csv")


## Inputs & parameters
inputs <- paste0("\n MR Inputs: \n -pheno ", pheno, "\n -tag ", tag,  
                 "\n -covarSet ", covarSet, "\n -outcome ", outcome, 
                 "\n -ssInput ", ssInput, "\n -datInput ", datInput, 
                 "\n -pctBthold ", pctBthold, "\n -tabBchangeAllCov ", tabBchangeAllCovInput,
                 "\n -mrDir ", mrDir)
cat(inputs)


# ==========================================
## Load data input files
# ==========================================

## Input files
ss<-fread(ssInput) %>% filter(LOCI == 1)
dat<-readRDS(datInput)

## tabBchangeAllCov
tabBchangeAllCov <- fread(tabBchangeAllCovInput)


## list of snps
snps <- names(dat %>% select(contains(ss$SNP)))

## covariate sets
covarSet <- covarSets[[covarSet]]

## outcome
outcome <- outcomes.l[[outcome]]


## Print covariate & outcomes 
cat("\n Running MR component of rediMR for\n pheno:", pheno, "\n on outcome:", outcome$Label,  " \n adjusting for covariate set:", 
    covarSet$Label, "\n", 
    paste0(" - ", t(covarSet$Names), "\n"))


##################################################
##  STEP 2: Two-Sample Mendelian Randomization  ##
##################################################

# install twosampleMR packages
#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github('MRCIEU/ieugwasr')
library(TwoSampleMR) ; library(ieugwasr)


# establish token for ieugwas access
opengwas_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqdWxpZS5nZXJ2aXNAdHVmdHMuZWR1IiwiaWF0IjoxNzE4NzM3MDA2LCJleHAiOjE3MTk5NDY2MDZ9.HLILykKTdeZJFp9-cYrNYtXR2xaTzQqPux5ojo3FTVHfWnJgkaSPozmeLHRYerD1sn8TSzc7l_FnPF4vFLPz7nufnp5tEtk4XButsDl9YXnjW-WvFgJ_sDDM1JucDi6IiO0-Bn9o5bppto9clUGku67sp3XYbwxxhQnCBCv-mfos8XI4x42y4GdjePOOklfLUiHxT_mt6F6tHKApyj0ssLC8T4b6cbEfD30jyA02V5XEhuUjSgFKM83rY6Qv0aB_SeaXygAV6sf5MkX4PTO49y7pcD4lDODqBwZlAroSWLikcWk1qKmfQvEQx5vadoyzHCcZlC5ATrz1Y5EMppAhhA"

cat("\n Starting Step 2: Two-Sample Mendelian Randomization ... \n ")


# ======================================================
## Build exposure_dat file for TwoSampleMR
# ======================================================

# -tabBchangeAllCov for refined/included SNP subsets
# -loci.afreq for effect allele frequency
# -ssInput.csv for [unrefined] gwas summary stats
afreq <- fread(paste0("../data/processed/rediMR/", pheno_tag, "_loci.afreq")) %>%
  mutate(REF_FREQS = 1-ALT_FREQS) %>% 
  select(SNP=ID, ALT_FREQS, REF_FREQS, REF_FREQ_ALLELE=REF, ALT_FREQ_ALLELE=ALT) 

exposure_dat <- ss %>% rowwise() %>% 
  mutate(refined.set = tabBchangeAllCov$RefinedSet[which(startsWith(tabBchangeAllCov$snp, SNP))]) %>%
  mutate(BETA_adj = tabBchangeAllCov$B_adj[which(startsWith(tabBchangeAllCov$snp, SNP))])  %>%
  mutate(SE_adj = tabBchangeAllCov$SE_adj[which(startsWith(tabBchangeAllCov$snp, SNP))])  %>%
  ungroup() %>% 
  mutate(SNP=gsub("snp", "", gsub("[.]", ":", SNP))) %>%
  left_join(afreq, by = "SNP") %>% 
  mutate(exposure = paste0(pheno), id.exposure=paste0("ukb-", pheno),
         effect_allele.exposure = ifelse(BETA>0, EA, NEA),
         other_allele.exposure = ifelse(BETA>0, NEA, EA),
         beta.exposure=abs(BETA), se.exposure=SE,
         beta.adj.exposure=abs(BETA_adj), se.adj.exposure=SE_adj) %>%
  mutate(eaf.exposure=ifelse(effect_allele.exposure==ALT_FREQ_ALLELE, ALT_FREQS, REF_FREQS)) %>%
  select(SNP, exposure, id.exposure, beta.exposure, se.exposure, beta.adj.exposure, se.adj.exposure,
         effect_allele.exposure, other_allele.exposure, eaf.exposure, refined.set)

# compile as list by snp refinement set
exposure_dat.l <- lapply(list(c(1,0),1,0), function(set) exposure_dat %>% filter(refined.set %in% set))
snp_sets <- c("All", "Refined", "Unrefined")
exposure_dat.l <- lapply(1:3, function(set) exposure_dat.l[[set]] %>% mutate(snp_set=snp_sets[set]))



# ======================================
## Compile outcomes from ieuGWAS
# ======================================

#outcome.Labs <- c(sapply(1:length(outcomes.l), function(i) outcomes.l[[i]]$label))
#outcome.vars <- c(sapply(1:length(outcomes.l), function(i) outcomes.l[[i]]$outcome_var))

outcome_dat.l <- lapply(1:length(outcomes.l), function(i) {
  extract_outcome_data(snps=exposure_dat$SNP, outcome = outcomes.l[[i]]$id.outcome, rsq=0.5, opengwas_jwt = opengwas_jwt) %>%
    mutate(type = rep(outcomes.l[[i]]$type))
}) ; names(outcome_dat.l) <- outcome.Labs

nOutcomes=length(outcome_dat.l)


##################
##    Run MR    ##
##################

# Write function to run MR for an exposure_dat & outcome_dat
## Run over list of SNPs ------------------------------------------------

run_mr <- function(exposure_dat, outcome_dat, adj_beta) {
  
  if(adj_beta==T) {
    exposure_dat
  }
  
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


## If using adjusted betas, 
if(adjBetaMR == "adjB") {
  MRsum_file_pref <- paste0(pheno_tag_covarSet, "_MRsum_adjB_")
} else { MRsum_file_pref <- paste0(pheno_tag_covarSet, "_MRsum_") }

lapply(as.list(1:nOutcomes), function(i) {
  saveRDS(mr_results.l[[i]], file = paste0(mrDir, "/", MRsum_file_pref, names(outcome_dat.l)[i], ".rda"))
})



###########################################################
##  Run MR Sensitivity Analysis - Validity vs. Strength  ##
###########################################################




##########################
## STEP 2 COMPLETED !!  ##
##########################

cat(paste0("
Congratulatulations!!\n",
           "You completed ReDiMR Step 2: Two-Sample MR for ", pheno, " with covariate set ",
           adjCovarSets[[covarSet]]$adjCovarLabel, " on the following outcomes:"))
cat(c(sapply(outcomes, function(x) x$outcome_var)))

cat("Data files are now ready for visualization! Have Fun!")

################################################################################




#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#
#////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////#



#EOF