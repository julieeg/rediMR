# ReDiMR
# Last updated: May 17, 2023




###########################################
##  STEP 1: Set up & assign parameters  ##
##########################################

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
pheno_tag <- paste0(pheno, "_", tag)

ssInput <- args[3] #paste0("../data/processed/rediMR/", pheno_tag, "_ssInput.csv")
datInput <- args[4] #paste0("../data/processed/rediMR/", pheno_tag, "_datInput.rda") 
pctBthold <- args[5] #20 
outDir <- args[6] #dirname(datInput)


# load basic functions
source("../scripts/basic_functions.R")


# =======================
##  Assign parameters
# =======================

# covariates in base gwas
gwasCovars <- c("age","sex", paste0("gPC", 1:10))


# covariates to adjust for in ReDiMR
adjCovars <- c("smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab", 
               "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10))

## Format covariates for adjustment 
adjCovarNames <- c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    pa_met_excess_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist-to-hip",
    dietPC1="Diet Pattern PC1", dietPC2="Diet Pattern PC2", dietPC3="Diet Pattern PC3",
    dietPC4="Diet Pattern PC4", dietPC5="Diet Pattern PC5", 
    #`dietPC1+dietPC2+dietPC3+dietPC4+diePC5`="Top 5 Diet Patterns", 
    dietPC6="Diet Pattern PC6", dietPC7="Diet Pattern PC7", dietPC8="Diet Pattern PC8", 
    dietPC9="Diet Pattern PC9", dietPC10="Diet Pattern PC10"
)

# Write as lists
adjCovars.l <- as.list(adjCovars)
adjCovarNames.l <- as.list(adjCovarNames)

nCovars.l <- as.list(1:length(adjCovarNames.l))


## Inputs & parameters
cat("\n ReDiMR Inputs: \n -ssInput ", ssInput, "\n -datInput ", datInput, "\n -pctBthold ", pctBthold, "\n -outDir ", outDir, "\n")

cat("\n Covariates for adjustment: \n", paste0(" - ", t(list2DF(adjCovarNames.l)), "\n"))




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

# list of snps
snps <- names(dat %>% select(contains(ss$SNP)))


# =============================================
## Write functions to calculate %B change 
# =============================================

pctBchange.fun <- function(pheno, snp, adjCovar, baseCovars=gwasCovars, replace_covar_name=NA, data=dat) {
  
  covarsFormatted <- paste0(baseCovars, collapse="+")
  snp_EA <- names(data %>% select(starts_with(snp)))
  baseM <- lm(formula(paste0(pheno, "~", snp_EA, "+", covarsFormatted)), data)
  adjM <- lm(formula(paste0(pheno, "~", snp_EA, "+", covarsFormatted, "+", adjCovar)), data)
  
  baseCI <- confint(baseM)[2,]
  adjCI <- confint(adjM)[2,]
  pctB <-round(((adjM$coef[2]-baseM$coef[2])/baseM$coef[2])*100, 2)
  
  baseP=summary(baseM)$coef[2,4]
  adjP=summary(adjM)$coef[2,4]

  if(!is.na(replace_covar_name)) {covarName = replace_covar_name} else {covarName = adjCovar}
  lmOut <- cbind.data.frame(
    snp=snp_EA, covar=covarName, B_base=baseM$coef[2], B_adj=adjM$coef[2], 
    lowCI_base=baseCI[1], upCI_base=baseCI[2], lowCI_adj=adjCI[1], upCI_adj=adjCI[2], 
    P_base=baseP, P_adj=adjP, B_pctChange=pctB) ; rownames(lmOut) <- paste0(covarName)
  
  return(lmOut)
}




########################################################################
## Basic descriptives of FOOD & covariates for adjustment
########################################################################

print_summary_table <- function(pheno) {
  do.call(rbind.data.frame, lapply(1:length(adjCovarNames), function(i) {
    covar=adjCovarNames[i]
    tmp <- dat %>% select(Covariate=all_of(names(covar)), y=all_of(pheno)) %>% filter(!is.na(y)) 
    if(is.factor(tmp$Covariate)) {
      out <- tmp %>% group_by(Covariate) %>% 
        summarise(Mean.SD = mean_sd(y)) %>% mutate(Variable.Name=covar,.before=Covariate) %>% 
        as.data.frame() } 
    else if(is.numeric(tmp$Covariate)) {
      out <- tmp %>% mutate(Covariate=ifelse(Covariate<median(Covariate, na.rm=T), "Above median", "Below median")) %>%
        group_by(Covariate) %>% summarise(Mean.SD=mean_sd(y)) %>% 
        mutate(Variable.Name=covar, .before=Covariate) %>% as.data.frame() }
  } ))
}

print_summary_table("raw_veg") %>% write.csv(file = paste0(outDir, "/", pheno_tag, "_tabDescrByCov.csv"), row.names=T)




########################################################################
## Run rediMR to calculate pct Beta change after covariate adjustment ## 
########################################################################

# ================================================
## Adjust for ALL covariates --> SNP Refinement  
# ================================================

cat("\n Calculating pctBchange when adjusting for ALL covariates ... \n ")

# For eacn SNP, tabulate pctBchange when adjusting for ALL covariates
tabBchangeAllCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) {
  pctBchange.fun(pheno=pheno, snp=snp, adjCovar=paste0(adjCovars, collapse="+"), 
                 replace_covar_name = "All_Covariates", data=dat)}, mc.cores = 8 )) %>% 
  mutate(RefinedSet = ifelse(abs(B_pctChange) < as.numeric(pctBthold),1,0)) 


# Write results to csv
write.csv(tabBchangeAllCov, file = paste0(outDir, "/", pheno_tag, "_tabBchangeAllCov.csv"), row.names=T)

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno_tag, "_tabBchangeAllCov.csv")))
head(tabBchangeAllCov)



# =======================================================
## Adjusting for EACH covariate --> Covariate Influence  
# =======================================================

cat("Calculating pctBchange when adjusting for EACH covariate ... ")

# For each snp, tabulate pctBchange when adjusting for EACH covariate
tabBchangeByCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) { 
  do.call(rbind.data.frame, mclapply(nCovars.l, function(i) {
    pctBchange.fun(pheno=pheno, snp=snp, adjCovar=names(adjCovarNames)[i], 
                   replace_covar_name=adjCovarNames[[i]], data=dat) }, mc.cores = 8))  }, mc.cores = 8))

# Write results to csv
fwrite(tabBchangeByCov, file = paste0(outDir, "/", pheno_tag, "_tabBchangeByCov.csv"))

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno_tag, "_tabBchangeByCov.csv")))
head(tabBchangeByCov)


###############################################################################
                      ## ~~~ Step 1 Completed !! ~~~~ ##
cat(paste0("\n Congratulatulations!!\n",
          " You completed ReDiMR Step 1: SNP Refinement for ", pheno, "!\n"),
    "Now proceeding to Step 2: Two-Sample MR with", mr_traits, "...")
###############################################################################





###########################################
##  STEP 2: Set up & assign parameters  ##
##########################################

install.packages("ieugwasr")
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github('MRCIEU/ieugwasr')
library(TwoSampleMR) ; library(ieugwasr)


# establish token for ieugwas access
opengwas_jwt <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqdWxpZS5nZXJ2aXNAdHVmdHMuZWR1IiwiaWF0IjoxNzE3MTc5MzkxLCJleHAiOjE3MTgzODg5OTF9.Uu_WDLV7LmYN9ZYoEnVr_pOEUu9ggPHLH3gdUNJ4CcyBSG4VQ2eB-prv8SzaCINn8-pEExBBSnrUx6lHS5324g9qvQWV52cVZuwBEXjlvqxKz99Lq6AIRsUWaTxAZ38JWIiGF6M2bElty9G0NE4WPIWi76OTE9Ad4iv5SQZvFWO-4qwxzxudcWnsNasvd6Y12HSaTnMjDenR2DU8guF_O_bFL3yRxOjCKZT6WKfplmRZXJM8HADamK7h-FUUp7RMavgvLeLWjUhVT1pLVL3BvHSvOiH49lRl9FBG_KAuO6nnewWZCVh2eAof9L2Ud2kec4o_PIF8VsoAl8bwe-aUvA"


## Build exposure_dat file for TwoSampleME; using
# -tabBchangeAllCov for refined/included SNP subsets
# -loci.afreq for effect allele frequency
# -ssInput.csv for [unrefined] gwas summary stats

afreq <- fread(paste0("../data/processed/rediMR/", pheno_tag, "_loci.afreq")) %>%
  mutate(SNP=ifelse(!startsWith(ID, "rs"), gsub(":", ".", paste0("snp", ID)), ID)) %>%
  mutate(REF_FREQS = 1-ALT_FREQS) %>% 
  select(SNP, ALT_FREQS, REF_FREQS, REF_FREQ_ALLELE=REF, ALT_FREQ_ALLELE=ALT) 

exposure_dat <- ss %>% left_join(afreq, by = "SNP") %>% 
  mutate(effect.allele_exposure = ifelse(BETA>0, EA, NEA),
         other.allele_exposure = ifelse(BETA>0, NEA, EA),
         beta.exposure=abs(BETA), se.exposure=SE) %>%
  mutate(eaf.exposure=ifelse(effect.allele_exposure==ALT_FREQ_ALLELE, ALT_FREQS, REF_FREQS)) %>%
  rowwise() %>% mutate(refined.set = tabBchangeAllCov$RefinedSet[which(startsWith(tabBchangeAllCov$snp, SNP))]) %>%
  ungroup() %>% select(SNP, beta.exposure, se.exposure, effect.allele_exposure, other.allele_exposure, eaf.exposure, refined.set)



outcome_test=extract_outcome_data(snps=ss$SNP, outcomes = "ebi-a-GCST90018808", opengwas_jwt = opengwas_jwt)












