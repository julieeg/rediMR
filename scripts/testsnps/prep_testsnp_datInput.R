#load packages
library(tidyverse) ; library(data.table)



# command args
args <- commandArgs(trailingOnly = T)
pheno <- args[1]
tag <- "testsnps"
pheno_tag=paste0(pheno,"_",tag)

ANC="EUR"


testDir="../data/processed/testsnps"
phenoFile="../data/processed/ukb_phenos_unrelated.rda"



source("../scripts/basic_functions.R")




################################
## Relabel SNPs without rsIDs ##
################################

g<-fread(paste0(testDir, "/testsnps.raw"))

dat <- left_join(readRDS(phenoFile), g %>% rename(id=IID) %>% mutate(id=as.character(id)), by = "id") %>% filter(ancestry == ANC) 




#####################################
# Build diet PCs for diet patterns ##
#####################################

## Load in "raw" dietary data, without outlier removal (>5 SD)
diet_raw <- readRDS("../data/processed/ukb_phenos_raw_unrelated.rda")


if(pheno == "total_veg") { #gsub(".*_", "", pheno)
  vars_to_exclude <- c("raw_veg", "cooked_veg") 
} else if(pheno == "total_fruit") { 
  vars_to_exclude <- c("fresh_fruit", "dried_fruit") } else {
    vars_to_exclude <- paste0(pheno) 
  }


paste0("Building diet patterns via pca WITHOUT: ", vars_to_exclude)

vars_for_pca <- diet_raw %>% select(
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
  
  # Replace missing data with median values
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x)) %>%
  
  # Remove outliers
  mutate(across(where(is.numeric), function(i) remove_outliers.fun(i, SDs=5))) %>%
  filter(complete.cases(.)==T) %>%
  
  select(-(starts_with(vars_to_exclude)))


## Run PCA & save output as .rda
diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA
saveRDS(diet_pcs, paste0(testDir, "/", pheno_tag, "_dietPCs.rda"))


# Extract dietPC scores
dat.dietPCs <- cbind.data.frame(vars_for_pca$id, diet_pcs$x)
names(dat.dietPCs) <- c("id", paste0("dietPC", 1:ncol(diet_pcs$x)))


print("Done deriving diet patterns WITH:")
names(vars_for_pca)

dim(vars_for_pca)




################################################################################
## Merge data files & Exclude participants with >5SD for continuous variables ##
################################################################################

# Merge in diet data
dat.merged <- left_join(dat, dat.dietPCs, by = "id")

# Replace pheno values >5SD with NA **CONSIDER WINZORIZING**
#dat.merged <- dat.merged %>%
 # mutate(across(where(is.numeric), function(i) remove_outliers.fun(i, SDs=5)))



############################
## Merge data & write rda ##
############################

print("Compiling datasets and writing .rda file")

dat.merged %>% saveRDS(paste0(testDir, "/",pheno_tag,"_datInput.rda"))

print(paste0("The following variables are available for analysis: ", paste0(names(dat.merged), collapse=", ")))

##EOF



