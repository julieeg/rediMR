# ReDiMR - v1
# Last updated: May 15, 2023



############
## Set Up ##
############

# load required packages
library(tidyverse) ; library(data.table) ; library(dplyr) ; library(parallel)
library(paletteer) 

# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
ssInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno, "_ssInput.csv")  #args[2]
datInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno, "_datInput.rda") #args[3]
outDir <- dirname(datInput) #args[4]
pctBdeltIncl <- 20 #args[5]


# load basic functions
source("../scripts/basic_functions.R")

# create directory to store results
system(paste0("mkdir -p ", outDir))



## default ReDiMR variable assignments // optional arguments ?

# covariates in base gwas
gwasCovars <- c("age","sex", paste0("gPC", 1:10))

# covariates to adjust for in ReDiMR
adjCovars <- c("smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab", 
               "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10))



######################################
## Format covariates for adjustment ##
######################################

# Assign covar names
adjCovarNames <- c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    pa_met_excess_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist-to-hip",
    dietPC1="Diet Pattern PC1", dietPC2="Diet Pattern PC2", dietPC3="Diet Pattern PC3",
    dietPC4="Diet Pattern PC4", dietPC5="Diet Pattern PC5", dietPC6="Diet Pattern PC6", 
    dietPC7="Diet Pattern PC7", dietPC8="Diet Pattern PC8", dietPC9="Diet Pattern PC9", 
    dietPC10="Diet Pattern PC10"
)

# Write as lists
adjCovars.l <- as.list(adjCovars)
adjCovarNames.l <- as.list(adjCovarNames)

nCovars.l <- as.list(1:length(adjCovarNames.l))

############################################################
## Build function to calculate pctBchange from Base Model ##
############################################################

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
    snp=paste(snp), covar=covarName, B_base=baseM$coef[2], B_adj=adjM$coef[2], 
    lowCI_base=baseCI[1], upCI_base=baseCI[2], lowCI_adj=adjCI[1], upCI_adj=adjCI[2], 
    P_base=baseP, P_adj=adjP, B_pctChange=pctB) ; rownames(lmOut) <- paste0(covarName)
  
  return(lmOut)
}



###############################################
## Load input files & Run basic descriptives ##
###############################################

# Load summary stats data
ss<-fread(ssInput) %>% filter(LOCI == 1)

# Load phenotype/dosage data
dat<-readRDS(datInput)
head(dat)

# list of snps
snps <- names(dat %>% select(contains(ss$SNP)))


#######################################
## Basic descriptives by pheno_quant ##
#######################################

tab_descrByPhenoQ <- do.call(rbind.data.frame, lapply(nCovars.l, function(i) {
  tmp <- dat %>% select(phenoQ, cov=all_of(names(adjCovarNames)[i])) %>%
    filter(!is.na(phenoQ)) %>% group_by(phenoQ)
  if(is.factor(tmp$cov)) {
    out <- do.call(rbind.data.frame, lapply(levels(tmp$cov), function(lvl) {
      tmp %>% summarise(cov_level = n_pct(cov, level = lvl)) %>% 
        t() %>% as.data.frame() %>% filter(V1 != "Q1") } ))
    rownames(out) <- paste0( adjCovarNames[i], "_", levels(tmp$cov), ", n (%)") } 
  else if(is.numeric(tmp$cov)) {
    out <- tmp %>% summarise(cov_msd = mean_sd(cov)) %>%
      t() %>% as.data.frame() %>% filter(V1 != "Q1")
    rownames(out) <- paste0(adjCovarNames[i], ", mean \u00B1 SD") }
  colnames(out) <- c(paste0("Q", 1:ncol(out)))
  out
} )) ; head(tab_descrByPhenoQ)


cat(paste0("Writing file with descriptive characteristics by ", pheno, " quantiles: \n"))
print(dat %>% group_by(phenoQ) %>% select(phenoQ, pheno=all_of(pheno)) %>% summarise(m_SD=mean_sd(pheno, d=3)))

# Write descriptives to csv
fwrite(tab_descrByPhenoQ, file = paste0(outDir, "/", pheno, "_descrByPhenoQs.csv"))



############################################################
## Calculate pct % beta for ALL covariates --> Refinement ##
############################################################

cat("Calculating pctBchange when adjusting for ALL covariates ... ")

# For eacn SNP, tabulate pctBchange when adjusting for ALL covariates
tab_pctBchangeAllCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) {
  pctBchange.fun(pheno=pheno, snp=snp, adjCovar=paste0(adjCovars, collapse="+"), replace_covar_name = "All_Covariates", data=dat) 
}, mc.cores = 8 )) %>% mutate(RefinedSet = ifelse(abs(B_pctChange) < pctBdeltIncl,1,0)) 


# Write results to csv
fwrite(tab_pctBchangeAllCov, file = paste0(outDir, "/", pheno, "_pctBchangeAllCov.csv"))

cat (paste0("Done. Result written to ", paste0(outDir, "/", pheno, "_pctBchangeAllCov.csv")))
head(tab_pctBchangeAllCov)



###############################################################
## Calculate pct % beta BY covariate --> Covariate Influence ##
###############################################################

cat("Calculating pctBchange when adjusting for EACH covariate ... ")

# For each snp, tabulate pctBchange when adjusting for EACH covariate
tab_pctBchangeByCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) { 
  do.call(rbind.data.frame, mclapply(nCovars.l, function(i) {
    pctBchange.fun(pheno=pheno, snp=snp, adjCovar=names(adjCovarNames)[i], replace_covar_name=adjCovarNames[[i]], data=dat)
  }, mc.cores = 8))  }, mc.cores = 8))


# Write results to csv
fwrite(tab_pctBchangeByCov, file = paste0(outDir, "/", pheno, "_pctBchangeByCov.csv"))



