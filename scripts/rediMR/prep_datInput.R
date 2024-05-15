#load packages
library(tidyverse) ; library(data.table)


# command args
args <- commandArgs(trailingOnly = T)
pheno <- args[1]
datInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno, "_datInput.tmp") # args[2]
ssInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno, "_ssInput.csv") 
datOutput <- dirname(datInput)


## Load datInput.tmp file
dat <- fread(datInput)
ss <- fread(ssInput)

source("../scripts/basic_functions.R")


#########################
## Compile dosage data ##
#########################

dat.geno <- dat %>% select(id, contains(ss$SNP)) %>%
  # Rename snps within CHR:POS_REF_ALT coding to rsCHR.POS_REF_ALT
  rename_with(., ~ gsub(":", ".", ifelse(!startsWith(.x, "rs"), paste0("rs", .x), .x))) %>%
  rename(id=rsid)

ss %>% mutate(SNP = gsub(":", ".", ifelse(!startsWith(SNP, "rs"), paste0("rs", SNP ), SNP ))) %>%
  write.csv(ssInput, row.names=F)



############################################
##  recode character variables as factors ##
############################################

## Load & compile ukb phenotype data  -------------------------
paste("Recoding character variables as meaningful factors ...")

## basic phenotypes  ---------------
covars <- c(pheno, "age", "sex", paste0("gPC", 1:10), "smoke_level.lab", "alch_freq.lab", 
            "pa_met_excess_level.lab", "income_level.lab", "educ_level.lab", "bmi", "waist2hip") 

dat.covars <- dat %>%   
  # Recode no answer/do not know/missing as missing
  mutate(
    sex = case_when(sex == 1 ~ "Male", sex == 0 ~ "Female"),
    smoke_level.lab = case_when(
      smoke.lab == "No answer" ~ as.character(NA),
      smoke.lab != "No answer" ~ as.character(smoke.lab),
      TRUE ~ as.character(NA)),
    alch_freq.lab = case_when(
      alch_freq.lab == "Prefer not to answer" ~ as.character(NA),
      alch_freq.lab != "Prefer not to answer" ~ as.character(alch_freq.lab),
      TRUE ~ as.character(NA)),
    income_level.lab=case_when(
      income == "Prefer not to answer" ~ as.character(NA),
      income == "Do not know" ~ as.character(NA),
      income != "Prefer not to answer" & income != "Do not know" ~ as.character(income),
      TRUE ~ as.character(NA)),
    educ_level.lab = case_when(
      educ_level.lab == "Prefer not to answer" ~ as.character(NA),
      educ_level.lab != "Prefer not to answer" ~ as.character(educ_level.lab),
      TRUE ~ as.character(NA))) %>% 
  # Add descritptive labels & levels
  mutate(
    income_level.lab = factor(income_level.lab, 
                              levels = c("Less than 18,000", "18,000 to 30,999",
                                         "31,000 to 51,999", "52,000 to 100,000",  
                                         "Greater than 100,000"),
                              labels = c("lt_18000", "from_18000_to_30999",
                                         "from_31000_to_51999", "from_52000_to_100000",  
                                         "gt_100000")),
    #Edu levels & yrs based on: https://www.nature.com/articles/s41380-019-0596-9#MOESM1)
    educ_level.lab = factor(educ_level.lab, 
                            levels = c("College or university degree", # ~20yrs 
                                       "NVQ/HND or equivalent", # 2 of 3 years bachelor's degree ~19yrs
                                       "Other professional qualifications", # e.g., nursing degree, teaching degree ~ 15yrs
                                       "A/AS levels or equivalent", # 1 year bachelor's degree ~13yrs
                                       "O/GCSE levels or equivalent", # HS + Associates degree ~10yrs
                                       "CSEs or equivalent", # completed HS ~10yrs
                                       "None of the above")), #~7yrs
    smoke_level.lab = factor(smoke_level.lab, levels = c("Current", "Former", "Never")),
    alch_freq.lab = factor(alch_freq.lab, ordered = T,
                           levels = c("Daily or almost daily",
                                      "3-4 per week", "1-2 per week", "1-3 per month",
                                      "Special occasions only", "Never")),
    pa_met_excess_level.lab = factor(pa_met_excess_lvl, ordered = T,
                               levels = c("Low", "Moderate", "High"))) %>% 
  select(id, ancestry, c(paste0("gPC", 1:10)), all_of(covars))



#####################################
# Build diet PCs for diet patterns ##
#####################################

if(gsub(".*_", "", pheno) == "veg") { vars_to_exclude <- c("raw_veg", "cooked_veg") }
if(gsub(".*_", "", pheno) == "fruit") { vars_to_exclude <- c("fresh_fruit", "dried_fruit") }

paste0("Building diet patterns via PCA WITHOUT: ", paste0(vars_to_exclude, collapse = " and "))

vars_for_pca <- dat %>% select(
  id,
  cooked_veg_QT=cooked_veg, raw_veg_QT=raw_veg,
  fresh_fruit_QT=fresh_fruit, dried_fruit_QT=dried_fruit, 
  oily_fish_QT=oily_fish, nonoily_fish_QT=nonoily_fish,
  procmeat_QT=procmeat, poultry_QT=poultry, cheese_QT=cheese,
  beef_QT=beef, lamb_QT=lamb, pork_QT=pork, 
  bread_type_white_vs_brown_or_whole_BIN=bread_type_white_vs_brown_or_whole, 
  bread_intake_QT=bread_intake,
  milk_type_full_vs_low_or_nonfat_BIN=milk_type_full_vs_low_or_nonfat,
  cereal_type_sugar_vs_any_bran_BIN=cereal_type_sugar_vs_any_bran, cereal_intake_QT=cereal_intake,
  spread_type_butter_vs_any_other_BIN=spread_type_butter_vs_any_other,
  coffee_type_decaf_vs_regular_BIN=coffee_type_decaf_vs_regular, coffee_QT=coffee,
  tea_QT=tea, water_QT=water, 
  addsalt_always_often_vs_nrs_BIN=addsalt_always_often_vs_nrs,
  hotdrink_temp_hot_or_vhot_vs_warm_BIN=hotdrink_temp_hot_or_vhot_vs_warm) %>%
  
  # Median impute for outliers >5sd
  mutate(across(where(is.numeric), function(i) remove_outliers.fun(i, SDs=5))) %>%
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x))  %>%
  
  select(-starts_with(vars_to_exclude))


## Run PCA & save output as .rda
diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA
saveRDS(diet_pcs, paste0(datOutput, "/", pheno, "_dietPCs.rda"))


# Extract dietPC scores
dat.dietPCs <- cbind.data.frame(vars_for_pca$id, diet_pcs$x)
names(dat.dietPCs) <- c("id", paste0("dietPC", 1:ncol(diet_pcs$x)))


# Compile diet variables 
dat.diet <- left_join(vars_for_pca, dat.dietPCs, by = "id")


print("Done compiling diet variables.")



################################################################################
## Merge data files & Exclude participants with >5SD for continuous variables ##
################################################################################

dat.merged <- dat.covars %>% 
#  left_join(dat.phenoQ, by = "id") %>%
  left_join(dat.diet, by = "id") %>%
  left_join(dat.geno, by = "id") 


# Replace with NA values >5SD for pheno
dat.merged <- dat.merged %>%
  mutate(across(pheno, ~ remove_outliers.fun(.x, SDs=5))) 


#################################################
## Create quantile variable for diet phenotype ##
#################################################

phenoQs <- quantile(dat.merged %>% select(pheno), probs=seq(0,1,0.25), na.rm=T)
if(length(unique(phenoQs)) != length(phenoQs)) {
  phenoQs <- quantile(dat.merged %>% select(pheno), probs=seq(0,1,0.33), na.rm=T)
}

dat.phenoQ <- dat.merged %>% select(id, phenoX=all_of(pheno)) %>%
  mutate(phenoQ = cut(phenoX, breaks = phenoQs, labels = c(paste0("Q", 1:(length(phenoQs)-1))), include.lowest=T)) %>%
  mutate(phenoQ = factor(phenoQ, levels = c(paste0("Q", 1:(length(phenoQs)-1))))) %>%
  select(id, phenoQ) #%>% 


############################
## Merge data & write rda ##
############################

print("Compiling datasets and writing .rda file")

dat.merged <- dat.merged %>% 
  left_join(dat.phenoQ, by = "id") 

dat.merged %>% saveRDS(paste0(datOutput, "/", pheno, "_datInput.rda"))


print(paste0("DONE merging phenotype/dosage files for ", pheno))
print(paste0("The following variables are available for analysis: ", paste0(names(dat.merged), collapse=", ")))

##EOF



