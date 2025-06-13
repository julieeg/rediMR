## Prepare UKB phenotype data
## Last updated: 20205-06-04


############
## Set Up ##
############

# load basic packages & functions for data prep & cleaning
library(tidyverse) ; library(data.table)
source("../scripts/basic_functions.R")


# ========================
## Winsorize data by SD
# ========================

winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


# =========================
## Add descriptive labels 
# =========================

#e.g., labs_vals = female_labs <- list("Female" = 0, "Male" = 1)
descr_label.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}


descr_label_ordered.fun <- function(data, base_var, labs_vals) {
  base <- data %>% select(all_of(base_var)) 
  temp <- rep(NA, length(base))
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; temp <- factor(temp, levels=names(labs_vals)) 
  return(temp)
}


#######################################################
## Load & prepare demographic & Lifestyle variables  ##
#######################################################

print("Gathering demographic and lifestyle variables ...")


### Basic phenotypes ------------------

female_labs <- list("Female" = 0, "Male" = 1)
smoke_labs <- list("No answer"=-3, "Never"=0, "Former"=1, "Current"=2)
smoking_num <- c("0" = 0, "1" = 1, "2" = 2, "-9" = -3)
med_mets_labs <- list("Cholesterol lowering" = 1, "Blood pressure" = 2, "Insulin" = 3,
                      "None of the above" = -7, "Do not know" = -1, "Prefer not to answer" = -3)

base_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE)

base_phenos_id <- base_phenos %>% 
  select(id = f.eid,
         ac = f.54.0.0,
         ac_date = f.53.0.0,
         sex = f.31.0.0,
         age = f.21022.0.0,
         bmi = f.21001.0.0,
         waist = f.48.0.0,
         hip = f.49.0.0,
         smoking = f.20116.0.0) %>%
  mutate(
    female.lab = descr_label.fun(., "sex", female_labs),
    smoke.lab = descr_label.fun(., "smoking", smoke_labs),
    smoking.num = descr_label.fun(., "smoking", smoking_num),
    waist2hip = waist/hip) %>%
  mutate(
    smoke_level.lab = factor(case_when(
      smoke.lab == "No answer" ~ as.character(NA),
      smoke.lab != "No answer" ~ as.character(smoke.lab),
      TRUE ~ as.character(NA)),
      levels = c("Current", "Former", "Never"))
    )

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw/withdraw27892_232_14_Nov_2022.txt", what=character())


## Add Assessment Center ------------------

ac_labs <- list("Barts"=11012, "Birmingham" = 11021, "Bristol" =	11011, "Bury" =	11008, 
                "Cardiff" =	11003, "Cheadle (revisit)" =	11024, "Croydon" =	11020, 
                "Edinburgh" =	11005, "Glasgow" = 11004, "Hounslow" = 11018, "Leeds" = 11010,
                "Liverpool"=11016, "Manchester"=11001, "Middlesborough"=11017, "Newcastle" =11009, 
                "Nottingham"=11013, "Oxford"=11002, "Reading"=11007, "Sheffield"=11014, "Stockport (pilot)"=10003,
                "Stoke"=11006, "Swansea"=	11022,"Wrexham" =11023, "Cheadle (imaging)"=11025,
                "Reading (imaging)"=11026, "Newcastle (imaging)" =11027, "Bristol (imaging)"=11028)

base_phenos_id <- base_phenos_id %>% 
  mutate(ac.f = descr_label.fun(., "ac", ac_labs))


### Education level ------------------

## Coding based on: Ge T., et al. Cerebral Cortex 2019;29(8): 3471-3481.

educ_level_labs <- list(
  "None of the above" = -7, 
  "Prefer not to answer"= -3, 
  "College or university degree" = 1, 
  "A/AS levels or equivalent" = 2, 
  "O/GCSE levels or equivalent" = 3, 
  "CSEs or equivalent" = 4,
  "NVQ/HND or equivalent" = 5, 
  "Other professional qualifications" = 6)

educ_isced_level_labs <- list("Level 5" = 1, "Level 3" = 2, "Level 2" = 3, 
                              "Level 2" = 4, "Level 5" = 5, "Level 4" = 6, "Level 1" = -7)  # NA = -3 or missing

educ_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672670.tab.gz",
                 data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, educ_level = f.6138.0.0) %>%
  mutate(educ_level.lab = descr_label.fun(., "educ_level", educ_level_labs),
         educ_isced.lab = descr_label.fun(., "educ_level", educ_isced_level_labs)) %>% 
  mutate(
    educ_level.lab = case_when(
      educ_level.lab == "Prefer not to answer" ~ as.character(NA),
      educ_level.lab != "Prefer not to answer" ~ as.character(educ_level.lab),
      TRUE ~ as.character(NA))) %>% 
  mutate(
    #Edu levels & yrs based on: https://www.nature.com/articles/s41380-019-0596-9#MOESM1)
    educ_level.lab = factor(educ_level.lab, levels = c(
      "College or university degree", # ~20yrs 
      "NVQ/HND or equivalent", # 2 of 3 years bachelor's degree ~19yrs
      "Other professional qualifications", # e.g., nursing degree, teaching degree ~ 15yrs
      "A/AS levels or equivalent", # 1 year bachelor's degree ~13yrs
      "O/GCSE levels or equivalent", # HS + Associates degree ~10yrs
      "CSEs or equivalent", # completed HS ~10yrs
      "None of the above")) #~7yrs
  )


### Alcohol intake frequency -------------------------------------------

alch_freq_labs <- list("Prefer not to answer" = -3, "Daily or almost daily" = 1, 
                       "3-4 per week" = 2, "1-2 per week" = 3, "1-3 per month" = 4, "Special occasions only" = 5, "Never" = 6) 
alch_num <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "-9" = -3)
drinker_status_labs <- c("Never" = 0, "Previous" = 1, "Current" = 2)


## updated: 09-24-2024
alch_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669173.tab.gz",
                 data.table=FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid,
         alch_freq = f.1558.0.0,
         alch_drinker_status = f.20117.0.0,
         alch_redwine_wk = f.1568.0.0,
         alch_wht_chm_wk = f.1578.0.0,
         alch_beer_wk = f.1588.0.0,
         alch_spirit_wk = f.1598.0.0,
         alch_fortwin_wk = f.1608.0.0,
         alch_othr_wk = f.5364.0.0) %>%
  
  ## Get alcohol drinker status (never/former/current)
  mutate_at("alch_drinker_status", ~ifelse(is.na(.) | . == -3, NA, .)) %>%
  mutate(alch_drinker_status.lab = descr_label.fun(., "alch_drinker_status", drinker_status_labs)) %>%
  mutate(alch_drinker_status.lab = factor(descr_label.fun(., "alch_drinker_status", c("Never" = 0, "Previous" = 1, "Current" = 2)), levels=c("Never", "Previous", "Current"))) %>%
  
  ## Recode #drinks per week: do not know/no answer asmissing & winsorise
  mutate(across(ends_with("_wk"), ~ifelse(is.na(.)==T | . %in% c(-1,-3), 0, .))) %>%
  mutate(across(ends_with("_wk"), ~winsorize(.))) %>%
  
  ## Make variable for curent/nondrinkers
  mutate(alch_currdrinker = ifelse(is.na(alch_drinker_status) | alch_drinker_status == -3, NA, 
                                   ifelse(alch_drinker_status==2, 1, 0) )) %>%
  mutate(alch_neverdrinker = ifelse(!is.na(alch_drinker_status) & alch_drinker_status == 0, 1, 0)) %>%
  mutate(alch_gm_per_wk = alch_redwine_wk*16.8 + alch_wht_chm_wk*16.8 + alch_beer_wk*16 + 
           alch_spirit_wk*8 + alch_fortwin_wk*14.08 + alch_othr_wk*12) %>%
  mutate(alch_drinks_per_week = alch_gm_per_wk / 14 ) %>%
  
  # prepare alcohol frequency variable as ordered factor & numeric
  mutate(alch_freq.lab = descr_label.fun(., "alch_freq", alch_freq_labs),
         alch_freq.num = descr_label.fun(., "alch_freq", alch_num)) %>%
  mutate(alch_freq.lab = factor(alch_freq.lab, levels= c(names(alch_freq_labs)[-1]) )) %>%
  select(id, alch_freq.num, alch_freq.lab, alch_drinker_status, alch_drinker_status.lab, alch_drinks_per_week, alch_gm_per_wk) 
  
 
### Physical Activity -------------------------------------------

pa_fields <- c("walking_dur", "walking_frq", "moderate_dur", "moderate_frq",
               "vigorous_dur", "vigorous_frq")

pa_vars <- c("physact_met_excess", "physact_level")

pa_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb671173.tab.gz", 
               data.table=FALSE, stringsAsFactors=FALSE) %>% select(
                 id = f.eid,
                 iqpa_met = f.22040.0.0,
                 
                 # activity duration (mins) & frequency (days/week) 
                 walking_dur = f.874.0.0, walking_frq = f.864.0.0, #duration of walks (min) & days/week of walks + 10 min
                 moderate_dur = f.894.0.0, moderate_frq = f.884.0.0, #duration of moderate activity (min) & days/week of moderate activity
                 vigorous_dur = f.914.0.0, vigorous_frq = f.904.0.0, #duration of vigorous activity (min) & days/week of vigorous activity
                 pa_type = f.6164.0.0
               ) %>%  
  
  # if walking_frq = -2 (Unable to walk) --> Recode to 0
  mutate_at("walking_frq", ~ifelse(. == -2, 0, .)) %>%
  
  # if activity duration <10 min/day --> Recode to 0
  mutate(across(ends_with("dur"), ~ifelse(.<10,0,.))) %>%
  
  # replace Do not know (-1), Prefer not to answer (-3) or missing (NA) with median
  mutate(across(ends_with("dur") | ends_with("frq") , ~ ifelse(.<0 | is.na(.)==T, 0, .))) %>%
  
  # multiply minutes per day per activity by excess MET score, per activity type
  mutate(met_excess_walk = ((walking_dur * 2.3)/60)*walking_frq,
         met_excess_mod = ((moderate_dur * 3.0)/60)*moderate_frq,
         met_excess_vig = ((vigorous_dur * 7.0)/60)*vigorous_frq) %>%
  
  # calculate excess met-hr/wk
  mutate(physact_met_excess = met_excess_walk + met_excess_mod + met_excess_vig) %>%
  mutate(physact_met_excess = winsorize(physact_met_excess))

physact_met_excess_lvls <- quantile(pa_id$physact_met_excess, probs=seq(0,1,0.33), include.lowest=F)[-1]

pa_id <- pa_id %>% mutate(
  physact_level = case_when(
    physact_met_excess < physact_met_excess_lvls[1] ~ "1",
    physact_met_excess >= physact_met_excess_lvls[1] & 
      physact_met_excess < physact_met_excess_lvls[2] ~ "2",
    physact_met_excess >= physact_met_excess_lvls[2] ~ "3")) %>%
  mutate(physact_level.lab = factor(physact_level, levels=c("1", "2", "3"), labels=c("Low", "Moderate", "High"))) %>%
  
  select(id, physact_met_excess, physact_level, physact_level.lab, iqpa_met)


### Income level ---------------------------------------------------

income_fields <- c(income = 738)
income_coding <- c(
  "1" = "Less than 18,000",
  "2" = "18,000 to 30,999",
  "3" = "31,000 to 51,999",
  "4" = "52,000 to 100,000",
  "5" = "Greater than 100,000",
  "-1" = "Do not know",
  "-3" = "Prefer not to answer"
)

income_id <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                   data.table = FALSE, stringsAsFactors = FALSE) 
income_id <- income_id %>%
  mutate(income = income_coding[as.character(f.738.0.0)]) %>%
  select(id=f.eid, income) %>%
  mutate(income_level.lab = case_when(income == "Do not know" ~ mode(NA),
                                  income == "Prefer not to answer" ~ as.character(NA),
                                  income != "Do not know" & income != "Prefer not to answer" ~ as.character(income))) %>%
  mutate(income_level.lab = factor(income_level.lab, levels = c(
    "Less than 18,000", "18,000 to 30,999", "31,000 to 51,999", "52,000 to 100,000", "Greater than 100,000"),
    labels = c("lt_18000", "from_18000_to_30999", "from_31000_to_51999", "from_52000_to_100000",  "gt_100000"))
  )
    


#################################################
##   COMBINE all base + covariate variables    ##
#################################################

phenos_id <- base_phenos_id %>% 
  full_join(educ_id, by = "id") %>%
  full_join(alch_id, by = "id") %>%
  full_join(pa_id, by = "id") %>%
  full_join(income_id, by = "id")


print(paste0("DONE: Basic phenotypes prepared for", nrow(base_phenos_id), " participants"))


# ====================================================================
## Prepare FFQ dietary data 
# ====================================================================

## Vectors of FFQ variables --------------------

intake_fields <- c(
  cooked_veg = 1289, raw_veg = 1299, fresh_fruit = 1309, dried_fruit = 1319,
  bread_intake = 1438, bread_type = 1448, water = 1528, milk_type = 1418, 
  spread_type = 1428, spread_type_nonbutter=2654, cereal_intake = 1458, cereal_type = 1468,
  addsalt=1478, tea=1488, coffee = 1498, coffee_type = 1508, hotdrink_temp = 1518
) 

freq_fields <- c(
  oily_fish = 1329, nonoily_fish = 1339, procmeat = 1349, poultry = 1359, cheese = 1408,
  beef = 1369, lamb = 1379, pork = 1389
)

ffq_fields<-c(intake_fields, freq_fields)
ffq_vars <- setNames(paste0("f.", ffq_fields, ".0.0"), names(ffq_fields))


## Function to convert frequency values to servings/day (from KEW) --------------------

ffq_freq_to_sev <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}


# function to recode negative values as meaninginful --------------------

neg_to_num <- function(x) {
  #x <- as.double(x) #"double" required to add values with decimals (previously, integer)
  case_when(
    x >= 0 ~ as.numeric(x), # -1 = "Do not know" ; -3 = "Prefer not to answer"
    x == -10 ~ 0.5, # -10 = "Less than 1 serving/day"
    TRUE ~ as.numeric(NA)
  )
}


# compile ffq variables & add total variables  --------------------

ffq_id <- base_phenos %>% select(id=f.eid, all_of(ffq_vars)) %>%
  mutate(whole_bread = case_when(
    bread_type == 3 ~ 1,
    bread_type %in% c(1, 2, 4) ~ 0,
    TRUE ~ as.numeric(NA)
  )) %>%
  mutate(across(names(freq_fields), ffq_freq_to_sev)) %>%
  mutate(across(names(intake_fields), neg_to_num)) %>%
  mutate(bread_intake = bread_intake / 7, # bread intake was provided in slices/week
         cereal_intake = cereal_intake / 7 # cereal intake was provided in bowls/week)
         
  ) %>% # Add FFQ vars for PCA
  
  mutate(
    bread_type_white_vs_brown_or_whole = case_when(      
      bread_type == 1 ~ 1, bread_type == 2 | bread_type == 3 | 
        bread_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_full_vs_low_or_nonfat = case_when(       
      milk_type == 1 ~ 1, milk_type == 2 | milk_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    milk_type_rare_never_BIN = case_when(
      milk_type == 6 ~ 1, milk_type != 6 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_butter_vs_any_other = case_when(          
      spread_type == 1 ~ 1, spread_type == 2 | spread_type == 3 ~ 0, TRUE ~ as.numeric(NA)),
    spread_type_rare_never_BIN = case_when(
      spread_type == 0 ~ 1, spread_type != 0 ~ 0, TRUE ~ as.numeric(NA)),
    cereal_type_sugar_vs_any_bran = case_when(
      cereal_type == 5 ~ 1, cereal_type != 5 ~ 0, TRUE ~ as.numeric(NA)),
    coffee_type_decaf_vs_regular = case_when(              
      coffee_type == 1 ~ 1, coffee_type == 2 | coffee_type == 3 | 
        coffee_type == 4 ~ 0, TRUE ~ as.numeric(NA)),
    addsalt_freq_QT = addsalt,        # INCLUDE
    addsalt_always_often_vs_nrs = case_when(
      addsalt == 3 | addsalt == 4 ~ 1, addsalt == 1 | addsalt == 2 ~ 0, TRUE ~ as.numeric(NA)),
    hotdrink_temp_hot_or_vhot_vs_warm = case_when(
      hotdrink_temp == 1  ~ 1, hotdrink_temp == 3 | hotdrink_temp == 2 ~ 0, TRUE ~ as.numeric(NA))
  )

ffq_id %>% saveRDS("../data/processed/ukb_ffq_processed.rda")



# ====================================================================
## Prepare 24HR dietary data 
# ====================================================================

## From Cole et al., Heritability:
#A detailed 24HR questionnaire in which a subset of participants answered over 200 questions 
#on specific foods and beverages consumed (with quantities) in the preceding 24-hour day. 
#The 24HR was implemented as a questionnaire for the final 70 K in-person baseline assessment
#center participants from 2009–2010 and emailed four times to 320 K participants who consented 
#to re-contact via email between February 2011 and April 2012. Approximately 200 K individuals 
#have at least one and up to five recorded 24HR questionnaires.

#Each questionnaire was filtered for credible estimates of total energy intake 
#[≥1,000 kJ (UKB field 100002) and ≤20 MJ for males and ≤18 MJ for females (UKB field 100026)],
#typical dietary intake (UKB fields 100020 and 20085), completion duration greater than or 
#equal to 5 min (UKB field 20082), and overall completion (UKB field 20081). Additionally, 
#the participant could not be pregnant within 1 year of taking the 24HR nor have a cancer 
#diagnosis within the previous year (UKB fields 3,140 and 40005). All 24HR questions were 
#converted into 1/0 for yes/no to consumption; each categorical response was coded similarly
#[e.g., UKB field 20086 for special diet was converted into six binary variables for each response
#(gluten-free, lactose-free, low calorie, vegetarian, vegan, and a combined vegetarian or vegan field)]. 
#24HR questions pertaining to quantity consumed were also included as continuous variables.


# ===================================
##  Load 24HR & helper data files 
# ===================================

cat("Loading 24HR data ...\n")

path_to_24hr="/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_jul_2024/ukbb_app27892_diet_07242024.csv"

# 24hr data
diet24hr=fread(path_to_24hr) %>% rename(id = eid)

## ids for participants who withdrew consent
withdraw=scan("/humgen/florezlab/UKBB_app27892/withdraw/w27892_20241217.csv", what=character()) #488

diet24hr <- diet24hr %>% filter(!id %in% withdraw)

codebook=readxl::read_excel("../run/ukb_24hr_codebook.xlsx")



##########################
##  24HR data cleaning  ##
##########################

## Used pre-calculated food group weights (g) by UKB, from 24-hour recalls, as described
# in https://link.springer.com/article/10.1007/s00394-021-02535-x (N=93 food groups)

## Data cleaning
# -restrict plausible total energy intakes (Cole: 1000 kj/day & <20,000kj [M] or <18000 kj [F])
# -restrict to typical_diet_yesterday
# -restrict to >= 5 min completion time
# -restrict to overall completion (yes)

# ======================================================================
## Restrictions: typical diet, plausible total energy, complete data
# ======================================================================

diet_valid_mean <- function(field, df) {
  # For a given nutrient:
  # - Select all columns for that nutrient
  # - Tabule number of instances (not missing)
  # - if >=2 instances, calculate mean; if not, code as NA
  diet_fields_df <- lapply(0:4, function(i) {
    field_id <- paste0("p", field,"_i", i) # diet field 
    valid_field <- paste0("p100026_i", i) # ukb field for valid 24hr (energy & sex)
    typical_field <- paste0("p100020_i", i) # ukb field for typical diet yesterday
    reason_field <- paste0("p20085_i", i) # ukb field for reason for atypical diet
    duration_field <- paste0("p20082_i", i) # ukb field for time duration of 24hr
    valid_24hr <- df %>% select(id, field=field_id, valid=valid_field, typical=typical_field, 
                                reason=reason_field, duration=duration_field) %>%
      mutate(field=ifelse(valid == "" & typical != "No" & reason == "" & !duration<5, field, NA)) %>% # replace with NA, if invalid 24hr 
      select(id, field_id=field) %>% rename_with(., ~gsub("field_id", field_id, .))
    valid_24hr
  }) %>% reduce(full_join, by = "id") %>%
    mutate(n_valid_24hr = rowSums(!is.na(across(paste0("p", field, "_i",0:4))))) %>%  # count n of valid 24hr
    mutate(field_mean=ifelse(n_valid_24hr>0, rowSums(across(paste0("p", field, "_i",0:4)), na.rm=T)/n_valid_24hr, NA)) %>%
    select(id, field_mean) %>%
    rename_with(., ~gsub("field_mean", paste0("p", field, "_mean"), .)) %>%
    pull(ends_with("mean"))
}


# Run over all 93 diet group fields
diet_fields <- c("26002", codebook %>% pull(Field.ID))
diet_means <- lapply(diet_fields, function(f) {
  diet24hr %>% mutate(
    field_mean=diet_valid_mean(f, .)) %>%
    rename_with(., ~gsub("field", paste0("p", f), .)) %>%
    select(id, ends_with("mean"))
}) %>% reduce(full_join, by = "id")

#diet_means_raw %>% fwrite("../data/processed/ukb_diet_mean24hr_raw.txt", sep="\t")


# ======================================================================
## Calculate intake of 13 food (g) and 2 beverage (mL) groups 
# ======================================================================

FoodGroups <- unique(codebook$FoodGroup.Name)

build_foodgroup <- function(food_group, df) {
  # For a given nutrient:
  # - Select all columns for that nutrient
  # - Tabule number of instances (not missing)
  # - if >=2 instances, calculate mean; if not, code as NA
  diet_vars_id <- codebook$Var.ID[codebook$FoodGroup.Name==food_group]
  diet_vars_id <- paste0(diet_vars_id, "_mean") # diet_id_mean 
  vars_dat <- df %>% select(id, all_of(diet_vars_id)) %>%
    mutate(fg_dat = rowSums(across(all_of(diet_vars_id)), na.rm=T)) %>%
    select(id, fg_dat) %>%
    rename_with(., ~gsub("fg_dat", paste0(food_group, "_mean"), .)) %>%
    pull(ends_with("mean"))
  return(vars_dat)
}

diet_means <- diet_means %>% mutate(
  fg_bevs_alcoholic=build_foodgroup("fg_bevs_alcoholic", .),
  fg_cereals=build_foodgroup("fg_cereals", .),
  fg_dairy_products=build_foodgroup("fg_dairy_products", .),
  
  fg_eggs=build_foodgroup("fg_eggs", .),
  fg_fats_spreads=build_foodgroup("fg_fats_spreads", .),
  fg_fish_dishes=build_foodgroup("fg_fish_dishes", .),
  
  fg_fruits=build_foodgroup("fg_fruits", .),
  fg_meat_products=build_foodgroup("fg_meat_products", .),
  fg_meat_substitutes=build_foodgroup("fg_meat_substitutes", .),
  
  fg_mixed_dishes=build_foodgroup("fg_mixed_dishes", .),
  fg_bevs_nonalch=build_foodgroup("fg_bevs_nonalch", .),
  fg_nuts_seeds=build_foodgroup("fg_nuts_seeds", .),
  
  fg_condiments=build_foodgroup("fg_condiments", .),
  fg_sweets_snacks=build_foodgroup("fg_sweets_snacks", .),
  fg_vegetables_potatoes=build_foodgroup("fg_vegetables_potatoes", .)
  
)


# ======================================================================
## Calculate intake of derived food/beverage groups for PCA analysis
# ======================================================================

pc_groups <- unique(codebook$DietPC.Group)

build_pca_groups <- function(pc_group, df) {
  # For a given group:
  # - Select all columns for group
  # - Tabulate number of instances (not missing)
  # - if >=2 instances, calculate mean; if not, code as NA
  diet_vars_id <- codebook$Var.ID[codebook$DietPC.Group==pc_group]
  diet_vars_id <- paste0(diet_vars_id, "_mean") # diet_id_mean 
  vars_dat <- df %>% select(id, all_of(diet_vars_id)) %>%
    mutate(fg_dat = rowSums(across(all_of(diet_vars_id)), na.rm=T)) %>%
    select(id, fg_dat) %>%
    rename_with(., ~gsub("fg_dat", paste0(pc_group, "_mean"), .)) %>%
    pull(ends_with("mean"))
  return(vars_dat)
}

diet_for_pca <- lapply(pc_groups, function(pc) {
  pc_mean <- paste0(pc, "_mean")
  pc.dat <- diet_means %>% mutate(pc=build_pca_groups(pc, .)) %>% select(id, pc) %>%
    rename_at("pc", ~gsub("pc", pc_mean, .))
}) %>% reduce(full_join, by="id")


diet_all <- full_join(diet_means, diet_for_pca, by="id")
diet_all %>% fwrite("../data/processed/ukb_diet24hr.csv")
diet_all %>% saveRDS("../data/processed/ukb_diet24hr.rda")


## Combine 24HR and FFQ diet datasets ------------------

diet_all_id <- diet_all %>% 
  left_join(., ffq_id, by = "id") %>%
  rename(nut_kcal=p26002_mean) %>%
  
  ## Replace values >5SD with NA********************
  mutate(across(c(starts_with("fg_") | starts_with("pc_") | "nut_kcal"), ~winsorize(., SDs=5))) %>%
  mutate(across(c(names(ffq_vars)), ~winsorize(., SDs=5)))



######################
## genetic ancestry ##
######################

### Compile Pan-UKBB genetic PCs (use to create European subset) --------------

anc_rel_id <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                    data.table=FALSE, stringsAsFactors=FALSE) %>% 
  mutate(f.eid = as.character(f.eid),
         unrelated = !related_return2442) %>%
  select(id=f.eid, ancestry=pop_return2442, unrelated,
         one_of(paste0("PC", 1:20, "_return2442"))) %>%
  rename_at(vars(contains("PC")), ~gsub("PC", "gPC", gsub("_return2442", "", .))) %>%
  mutate(id=as.integer(id))


print("Breakdown of available data by relatedness & ancestry from PanUKBB: ")
print(table(anc_rel_id$unrelated, anc_rel_id$ancestry))


######################################
##  Load in dietary traits from JC  ##
######################################

## dietary phenotypes from JC -------------------

diet_traits_fromJC <- fread("../data/processed/gwas/fromJC/BOLT_UKB_diet_genoQCEUR450K_phenotypes_ffq_QC_PCA_manuscript_12132018_florezconverstion.csv") %>%
  rename(id=Florez_FID, 
         oilyfish_QT=oilyfish_overallfreq.1329.average_QT,
         alch_glasspermonth_QT=anyalcohol_glassespermonth.derived.average_QT,
         bread_type_BIN=bread_typeused.1448.average_bin4)

cat("Adding:", names(diet_traits_fromJC %>% select(-"id")))



### Merge phenotypes -------------------------------------------

phenos <- phenos_id %>%
  left_join(diet_all_id, by="id") %>%
  left_join(anc_rel_id, by="id") %>%
  left_join(diet_traits_fromJC, by = "id") %>%
  filter(!(id %in% withdrawn_consent)) %>%
  mutate(id = format(id, scientific=FALSE)) %>%
  mutate(IID = id, ., before=id)



##############################
## write PROCESSED datasets ##
##############################

# Unrelated subsets
phenos %>%
  filter(unrelated == TRUE) %>% 
  fwrite("../data/processed/ukb_phenos_unrelated.csv", row.names = F, col.names = T)


# Unrealted & EUR subjects
phenos %>%
  filter(unrelated == TRUE) %>%
  filter(ancestry == "EUR") %>%
  saveRDS("../data/processed/ukb_phenos_unrelated_EUR.rda")


print("Done preparing UKB Phenotype data. Datasets are ready for analysis.")

## Print Data Dictionary
sink("../data/processed/phenos/ukb_phenos_str.txt")
str(phenos)
sink()


##EOF


