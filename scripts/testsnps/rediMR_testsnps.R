# "RediMR Simulation" for testsnps




##################################
##  Set up & assign parameters  ##
##################################

# =======================
##  Set up
# =======================

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno=args[1]
tag="testsnps"
pheno_tag=paste0(pheno,"_",tag)

ANC="EUR"
testDir="../data/processed/testsnps"


datInput=paste0(testDir, "/", pheno_tag, "_datInput.rda")


# load basic functions
source("../scripts/basic_functions.R")



# =======================
##  Assign parameters
# =======================

pctBthold=20

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





#################################################################
##  Load input files ; Write functions to calculate %B change  ##
#################################################################

# ==========================================
## Load data input files
# ==========================================

# Load phenotype/dosage data
dat<-readRDS(datInput)

# list of snps
snps <- names(dat %>% select(starts_with("rs")))



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
write.csv(tabBchangeAllCov, file = paste0(testDir, "/", pheno_tag, "_tabBchangeAllCov.csv"), row.names=T)

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
fwrite(tabBchangeByCov, file = paste0(testDir, "/", pheno_tag, "_tabBchangeByCov.csv"))

cat (paste0("DONE: Results written to ", paste0(testDir, "/", pheno_tag, "_tabBchangeByCov.csv")))
head(tabBchangeByCov)





## EOF

