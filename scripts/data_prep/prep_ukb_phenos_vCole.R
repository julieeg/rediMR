
############################
## Prepare phenotype data ##
############################

#load packages
library(tidyverse) ; library(data.table)

## Load ukb_phenos + diet_traits from JC
ukb_phenos <- fread("../data/processed/ukb_phenos_unrelated_EUR.txt") %>%
  select(-"V1")

## dietary phenotypes from JC
diet_traits_fromJC <- fread("../data/processed/gwas/fromJC/BOLT_UKB_diet_genoQCEUR450K_phenotypes_ffq_QC_PCA_manuscript_12132018_florezconverstion.csv") %>%
  rename(id=Florez_FID, 
         oilyfish_QT=oilyfish_overallfreq.1329.average_QT,
         alch_glasspermonth_QT=anyalcohol_glassespermonth.derived.average_QT,
         bread_type_BIN=bread_typeused.1448.average_bin4)
cat("Adding:", names(diet_traits_fromJC %>% select(-"id")))


##!## 
add_conf_id <- readRDS("../data/processed/ukb_conf_addn_09242024.rda")

## Merge & save as .txt file (for plink)
dat <- full_join(ukb_phenos %>% select(-c(pa_met_excess, pa_met_excess.lvl)), diet_traits_fromJC, by = "id") %>% 
  full_join(add_conf_id %>% select(-c(alch_freq.lab, alch_freq.num)), by = "id") 


## Add updated confounder variables --------------
dat <- dat %>%
  mutate(meds_bp = ifelse(
    sex == "Female" & meds_female.lab == "Blood pressure" | sex == "Male" & meds_male.lab == "Blood pressure", 1, 
    ifelse(sex == "Female" & !is.na(meds_female.lab) | sex == "Male" & !is.na(meds_male.lab), 0, NA))) %>%
  mutate(
    sbp_adj = ifelse(meds_bp == 1, sbp+15, sbp),
    dbp_adj = ifelse(meds_bp == 1, dbp+10, dbp)
  )
  #fwrite("../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits.txt", row.names = F)

# ===========================
## Prepare phenotype data
# ===========================

dat <- dat %>% 
  mutate(
    smoke_level.num = case_when(
      smoke_level.lab == "Current" ~ 3, 
      smoke_level.lab == "Former" ~ 2, 
      smoke_level.lab == "Never" ~ 1,
      TRUE ~ as.numeric(NA)),
    alch_freq.num = as.numeric(alch_freq.num),
    physact_level.num = case_when(
      physact_level.lab == "Low" ~ 1,
      physact_level.lab == "Moderate" ~ 2,
      physact_level.lab == "High" ~ 3,
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
    educ_isced.num = case_when(
      educ_isced.lab == "Level 1" ~ 1,
      educ_isced.lab == "Level 2"~ 2,
      educ_isced.lab == "Level 3"~ 3,
      educ_isced.lab == "Level 4" ~ 4,
      educ_isced.lab == "Level 5" ~ 5),
    educ_years.num = as.numeric(educ_years),
    alch_freq.num = case_when(
      alch_freq.lab == "Daily or almost daily" ~ 5,
      alch_freq.lab == "3-4 per week" ~ 4, 
      alch_freq.lab == "1-2 per week" ~ 3, 
      alch_freq.lab == "1-3 per month" ~ 2, 
      alch_freq.lab == "Special occasions only" ~ 1, 
      alch_freq.lab == "Never" ~ 0, 
      TRUE ~ as.numeric(NA))
  )



## Add Assessment Center ----------------------------------------------------------

descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}

ac_labs <- list("Barts"=11012, "Birmingham" = 11021, "Bristol" =	11011, "Bury" =	11008, 
                "Cardiff" =	11003, "Cheadle (revisit)" =	11024, "Croydon" =	11020, 
                "Edinburgh" =	11005, "Glasgow" = 11004, "Hounslow" = 11018, "Leeds" = 11010,
                "Liverpool"=11016, "Manchester"=11001, "Middlesborough"=11017, "Newcastle" =11009, 
                "Nottingham"=11013, "Oxford"=11002, "Reading"=11007, "Sheffield"=11014, "Stockport (pilot)"=10003,
                "Stoke"=11006, "Swansea"=	11022,"Wrexham" =11023, "Cheadle (imaging)"=11025,
                "Reading (imaging)"=11026, "Newcastle (imaging)" =11027, "Bristol (imaging)"=11028)

dat <- dat %>% mutate(ac.f = descr_label.fun(., "ac", ac_labs))

dat %>% fwrite("../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits_09292024.txt", row.names = F)



################################
## Prepare summary stats data ##
################################

fread("../data/processed/gwas/fromJC/BOLTlmm_UKB_genoQCEURN455146_v3_diet_oilyfish_overallfreq.1329.average_QT_BgenSnps_mac20_maf0.005_info0.6.gz") %>%
  rename(ID=SNP, ALT=ALLELE1, REF=ALLELE0, EAF=A1FREQ, P=P_BOLT_LMM) %>%
  fwrite("../data/processed/gwas/oilyfish_QT_vCole.gwas", row.names=F, col.names=T, sep=" ")


fread("../data/processed/gwas/fromJC/BOLTlmm_UKB_genoQCEURN455146_v3_diet_bread_typeused.1448.average_bin4_BgenSnps_mac20_maf0.005_info0.6.gz") %>%
  rename(ID=SNP, ALT=ALLELE1, REF=ALLELE0, EAF=A1FREQ, P=P_BOLT_LMM) %>%
  fwrite("../data/processed/gwas/bread_type_BIN_vCole.gwas", row.names=F, col.names=T, sep=" ")


