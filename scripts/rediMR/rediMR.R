# ReDiMR
# Last updated: May 17, 2023



# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# load pantry file with stored parameters, basic functions & ggplot templates
source("../scripts/pantry.R")



##################################
##  Set up & input parameters  ##
##################################

# ==========================
## Load command arguments
# ==========================

# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1] #diet phenotype
tag <- args[2] #version "tag"
pheno_tag <- paste0(pheno, "_", tag)

ssInput <- args[3] #paste0("../data/processed/rediMR/", pheno_tag, "_ssInput.csv")
datInput <- args[4] #paste0("../data/processed/rediMR/", pheno_tag, "_datInput.rda") 
covarSet <- args[5]
pctBthold <- args[6]
outDir <- args[7]

# Add second tag 
if(length(args) == 8) {
  tag2 = args[8]
}


# covariate sets
adjCovars.l <- as.list(adjCovarSets[[covarSet]]$adjCovars)
adjCovarNames.l <- as.list(adjCovarSets[[covarSet]]$adjCovarNames)
nCovars.l <- as.list(1:length(adjCovarNames.l))


## Inputs & parameters
cat("\n ReDiMR Inputs: \n -pheno ", pheno, "\n -tag ", tag,  "\n -ssInput ", ssInput, "\n -datInput ", datInput, "\n -covarSet ", covarSet, "\n -pctBthold ", pctBthold, "\n -outDir ", outDir, "\n")
cat("\n Covariates for adjustment: \n", paste0(" - ", t(list2DF(adjCovarNames.l)), "\n"))


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

print_summary_table("raw_veg") %>% write.csv(file = paste0(outDir, "/", pheno_tag, "_tabDescrByCov_", covarSet, ".csv"), row.names=T)



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
write.csv(tabBchangeAllCov, file = paste0(outDir, "/", pheno_tag, "_tabBchangeAllCov_", covarSet, ".csv"), row.names=T)

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno_tag, "_tabBchangeAllCov_", covarSet, ".csv")))
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
fwrite(tabBchangeByCov, file = paste0(outDir, "/", pheno_tag, "_tabBchangeByCov_", covarSet, ".csv"))

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno_tag, "_tabBchangeByCov_", covarSet, ".csv")))
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


## PRS input files

# All
snpset_ss %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag, "_", covarSet, "_all_prsInput"))

# Refined
snpset_ss %>% 
  filter(RefinedSet == 1) %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag, "_", covarSet, "_ref_prsInput"))

# Unrefined
snpset_ss %>% 
  filter(RefinedSet == 0) %>% 
  select(ID, A1=EA, BETA) %>% 
  write_tsv(file=paste0("../data/processed/prs/", pheno_tag, "_", covarSet, "_unref_prsInput"))



###############################################################################
                      ## ~~~ Step 1 Completed !! ~~~~ ##
cat(paste0("\n Congratulatulations!!\n",
          " You completed ReDiMR Step 1: SNP Refinement for ", pheno, "!\n"),
    "Now proceeding to Step 2: Two-Sample MR with", mr_traits, "...")
###############################################################################


#EOF


