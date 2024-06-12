#load packages
library(tidyverse) ; library(data.table)


# command args
args <- commandArgs(trailingOnly = T)
pheno <- args[1]
tag <- args[2]

pheno_tag <- paste0(pheno, "_", tag)


datInput <- paste0("../data/processed/rediMR/", pheno_tag, "_datInput.tmp.rda") # args[2]
ssInput <- paste0("../data/processed/rediMR/", pheno_tag, "_ssInput.tmp.csv") 
outDir <- dirname(datInput)


## Load Input files
dat <- readRDS(datInput)
ss <- fread(ssInput)

source("../scripts/basic_functions.R")


# =================================
## Relabel SNPs without rsIDs 
# =================================

dat <- dat %>% 
  rename_with( ~ifelse(!startsWith(.x, "rs"), gsub(":", ".", paste0("snp", .x, recycle0 = T)), .x), contains(ss$SNP))

ss %>% mutate(SNP = gsub(":", ".", ifelse(!startsWith(SNP, "rs"), paste0("snp", SNP ), SNP ))) %>%
  fwrite(gsub(".tmp", "", ssInput), row.names=F)



#####################################
# Build diet PCs for diet patterns ##
#####################################


if(pheno == "total_veg") { #gsub(".*_", "", pheno)
  vars_to_exclude <- c("raw_veg", "cooked_veg") 
  } else if(pheno == "total_fruit") { 
    vars_to_exclude <- c("fresh_fruit", "dried_fruit") } else {
      vars_to_exclude <- paste0(pheno) 
  }


paste0("Building diet patterns via pca WITHOUT: ", vars_to_exclude)


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
  
  # Replace missing data with median values
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x)) %>%
  
  # Winzorise
  mutate(across(where(is.numeric), function(i) winsorize(i, SDs=5))) %>%
  filter(complete.cases(.)==T) %>%
  
  select(-(starts_with(vars_to_exclude)))


## Run PCA & save output as .rda
diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA
saveRDS(diet_pcs, paste0(outDir, "/", pheno_tag, "_dietPCs.rda"))


# Extract dietPC scores
dat.dietPCs <- cbind.data.frame(vars_for_pca$id, diet_pcs$x)
names(dat.dietPCs) <- c("id", paste0("dietPC", 1:ncol(diet_pcs$x)))


print("Done deriving diet patterns WITH:")
c(names(vars_for_pca))


## Summarize missing and ranges:
print("Table of N missing values per diet trait: \n ")
cbind.data.frame(
    N_complete=sapply(vars_for_pca, function(col) sum(ifelse(is.na(col)==F,1,0)) ),
    N_missing=sapply(vars_for_pca, function(col) sum(ifelse(is.na(col)==T,1,0)) ) , 
        Range=sapply(vars_for_pca, function(col) range(col, na.rm=T)) %>% t()
  )

dim(vars_for_pca)




########################################################
##  Merge data & write rda files for rediMR pipeline  ##
########################################################

print("Compiling datasets and writing .rda file")


# Merge in diet data
dat.merged %>% 
  left_join(dat, dat.dietPCs, by = "id") %>%
  saveRDS(paste0(outDir, "/", pheno_tag, "_datInput.rda"))


print(paste0("DONE merging phenotype/dosage files for ", pheno))
print(paste0("The following variables are available for analysis: ", paste0(names(dat.merged), collapse=", ")))



##EOF



