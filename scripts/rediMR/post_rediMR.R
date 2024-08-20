#compared_rediMR

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
group <- args[3]


pheno_tag <- paste0(pheno, "_", tag)
pheno_tag_covarGroup <- paste0(pheno_tag, "_", group)

redimrDir <- "../data/processed/rediMR"
ssInput <- paste0(redimrDir, "/", pheno_tag, "_ssInput.csv")
datInput <- paste0(redimrDir, "/", pheno_tag, "_datInput.rda") 



## Inputs & parameters
inputs <- paste0("\n compare_rediMR Inputs: \n -pheno ", pheno, "\n -tag ", tag,  
                 "\n -group ", group)
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
covarGroups <- covarSetGroups[[group]]

cat(paste0("Comparing the following covariate groups:"), 
    paste0("\n - ", covarGroups$Sets), "\n")



################################################################################
## Confounder-phenotype associations, after adjusting for covariate sets
################################################################################

# =========================================================
# Add numeric versions of descriptive factor variables
# ==========================================================

dat <- dat %>% 
  mutate(
    smoke_level.num = case_when(
      smoke_level.lab == "Current" ~ 1, 
      smoke_level.lab == "Former" ~ 2, 
      smoke_level.lab == "Never" ~ 3,
      TRUE ~ as.numeric(NA)),
    alch_freq.num = as.numeric(alch_freq.num),
    pa_met_excess_level.num = case_when(
      pa_met_excess_level.lab == "Low" ~ 1,
      pa_met_excess_level.lab == "Moderate" ~ 2,
      pa_met_excess_level.lab == "High" ~ 3,
      TRUE ~ as.numeric(NA)),
    income_level.num = case_when(
      income_level.lab == "lt_18000" ~ 1,
      income_level.lab == "from_18000_to_30999" ~ 2,
      income_level.lab == "from_31000_to_51999" ~ 3,
      income_level.lab == "from_52000_to_100000" ~ 4,
      income_level.lab == "gt_100000" ~ 5,
      TRUE ~ as.numeric(NA)),
    educ_level.num = case_when(
      educ_level.lab == "None of the above" ~ 1,
      educ_level.lab == "CSEs or equivalent" ~ 2, # completed HS ~10yrs
      educ_level.lab == "O/GCSE levels or equivalent" ~ 3, # HS + Associates degree ~10yrs
      educ_level.lab == "A/AS levels or equivalent" ~ 4, # 1 year bachelor's degree ~13yrs
      educ_level.lab == "Other professional qualifications" ~ 5, # e.g., nursing degree, teaching degree ~ 15yrs
      educ_level.lab == "NVQ/HND or equivalent" ~ 6, # 2 of 3 years bachelor's degree ~19yrs
      educ_level.lab == "College or university degree" ~ 7, # 20yers+
      TRUE ~ as.numeric(NA)),
    educ_years.num = as.numeric(educ_years)
  )

head(dat %>% select(all_of(confounders.num)))



# ===========================================================================
## Examine confounder-pheno associations, after covariate adjustment
# ===========================================================================

## Unadjusted ============

#tab_CorConfdPheno <- do.call(
#  
#  rbind.data.frame, lapply(1:length(confounders.num), function(i) {
#    base=summary(lm(formula(paste0(pheno, "~", confounders.num[[i]])), data=dat))
#    
#  do.call(rbind.data.frame, lapply(as.list(covarGroups$Sets), function(set) {
#    covSet = covarSets[[set]]
#    adj=summary(lm(formula(paste0(pheno, "~", confounders.num[[i]], "+", paste0(covSet$Covars, collapse="+"))), data=dat))
#    out = cbind.data.frame(values=adj$coef[2,]) %>%
#      mutate(stats=c("B", "SE", "Tstat", "P")) %>% 
#      pivot_wider(values_from = "values", names_from=stats) %>% 
#      mutate(adjR2 = adj$adj.r.squared, Fstat=adj$fstat[1]) %>%
#      mutate(covar=covSet$Label, .before="B")
#    return(out) }) ) %>% 
#    bind_rows(cbind.data.frame(
#      covar="Unadjusted", B=base$coef[2,1], SE=base$coef[2,2], Tstat=base$coef[2,3],
#      P=base$coef[2,4], adjR2=base$adj.r.squared, Fstat=base$fstat[1])) %>%
#      mutate(
#        lowCI=B-1.96*SE, 
#        upCI=B+1.96*SE) %>%
#      mutate(confName=paste0(names(confounders.num)[[i]]), 
#             conf=paste0(confounders.num[[i]]), .before="covar")
#    } )
#  )

#tab_CorConfdPheno %>% saveRDS(paste0(redimrDir, "/", pheno_tag_covarGroup, "_tabCorConfdPheno.rda"))



## With sex- & age- adjustment ============

#base=summary(lm(formula(paste0(pheno, "~sex+age")), data=dat))

#tab_CorConfdPheno_sexage <- do.call(
#  rbind.data.frame, lapply(1:length(confounders.num), function(i) {
#    base=summary(lm(formula(paste0(pheno, "~", confounders.num[[i]], "+sex+age")), data=dat))
#    do.call(rbind.data.frame, lapply(covarGroups$Sets, function(set) {
#      covSet = covarSets[[set]]
#      adj=summary(lm(formula(
#        paste0(pheno, "~", confounders.num[[i]], "+sex+age+", paste0(covSet$Covars, collapse="+"))), data=dat))
#      out = cbind.data.frame(values=adj$coef[2,]) %>%
#        mutate(stats=c("B", "SE", "Tstat", "P")) %>% 
#        pivot_wider(values_from = "values", names_from=stats) %>% 
#        mutate(adjR2 = adj$adj.r.squared, Fstat=adj$fstat[1]) %>%
#        mutate(covar=covSet$Label, .before="B")
#      return(out) }) ) %>% 
#      bind_rows(cbind.data.frame(
#        covar="Unadjusted", B=base$coef[2,1], SE=base$coef[2,2], Tstat=base$coef[2,3],
#        P=base$coef[2,4], adjR2=base$adj.r.squared, Fstat=base$fstat[1])) %>%
#      mutate(
#        lowCI=B-1.96*SE, 
#        upCI=B+1.96*SE) %>%
#      mutate(confName=paste0(names(confounders.num))[[i]], 
#             conf=paste0(confounders.num[[i]]), .before="covar")
#  } )
#)

#tab_CorConfdPheno_sexage %>% saveRDS(paste0(redimrDir, "/", pheno_tag_covarGroup, "_tabConfdPheno_sexage.rda"))



################################################################################
## Confounder-phenotype associations, after adjusting for covariate sets
################################################################################


# ======================================================
## Heat Map of pheno-dietPC-confounder correlations
# ======================================================

# Create correlations matrix of all dietPCs with confounder

#cordat <- dat %>% select(pheno, starts_with(c("dietPC")), all_of(confounder_vars.num)) %>% filter(complete.cases(.))
#cormat <- cor(
#  cordat %>% select(pheno, starts_with("dietPC")), 
#  cordat %>% select(pheno, confounder_vars.num), method = "pearson")
#colnames(cormat) <- c(pheno, names(confounders.num))

#cormat %>% saveRDS(paste0("../data/processed/plots/", pheno_tag_covarGroup, "_plotHMconf.rda"))



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



