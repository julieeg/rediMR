## Prepare dietary phenotypes for diet PC exploratory analysis

library(tidyverse) ; library(data.table)


# ========================
## Winsorize data by SD
# ========================

winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


##################################################
##  Dietary data and diet PCs from FFQs & 24HR  ## 
##################################################

print("Preparing dietary data ...")

### Basic phenotypes ------------------

female_labs <- list("Female" = 0, "Male" = 1)

base_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE)


# ==================================
## Prepare FFQ dietary data 
# ==================================

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


# compile ffq variables & add total variables    
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


# ==================================
## Prepare 24HR dietary data 
# ==================================

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
diet_means_raw <- lapply(diet_fields, function(f) {
  diet24hr %>% mutate(
    field_mean=diet_valid_mean(f, .)) %>%
    rename_with(., ~gsub("field", paste0("p", f), .)) %>%
    select(id, ends_with("mean"))
}) %>% reduce(full_join, by = "id")

#diet_means_raw %>% fwrite("../data/processed/ukb_diet_mean24hr_raw.txt", sep="\t")

diet_means <- diet_means_raw %>%
  # Add total energy variable
  mutate(across(-id, ~winsorize(., SDs=5)))


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

build_pc_groups <- function(pc_group, df) {
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

food_groups_for_pca <- lapply(pc_groups, function(pc) {
  pc_mean <- paste0(pc, "_mean")
  pc.dat <- diet_means %>% mutate(pc=build_pc_groups(pc, .)) %>% select(id, pc) %>%
    rename_at("pc", ~gsub("pc", pc_mean, .))
  }) %>% reduce(full_join, by="id")


diet_all <- full_join(diet_means, food_groups_for_pca, by="id")
diet_all %>% fwrite("../data/processed/ukb_diet24hr.csv")
diet_all %>% saveRDS("../data/processed/ukb_diet24hr.rda")


#######################################
##  Prepare diet PCs from 24HR data  ##
#######################################

# Build diet PCs for diet patterns -------------

vars_for_pca <- diet_all %>% 
  mutate(id=as.character(id)) %>% 
  #select(id, all_of(FoodGroups)) %>%
  select(id, nut_kcal="p26002_mean", all_of(paste0(pc_groups, "_mean"))) %>%
  rename_with(., ~gsub("_mean", "_adjkcal", .)) %>%
  
  # Replace missing values with medians
  mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x)) %>%
  
  mutate_at(paste0(pc_groups, "_adjkcal"), ~ resid(lm(.x ~ nut_kcal), data=.)) %>%
  
  # Remove nut_kcal
  select(-nut_kcal)
  
  # Winsorize data to 5 SD
  #mutate(across(where(is.numeric), function(i) winsorize(i, SDs=5))) %>%

## Run PCA
set.seed(314159)
dietPCs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA

# compile dietPC scores
dietPCs_id <- as.data.frame(cbind(id=vars_for_pca$id, dietPCs$x))
colnames(dietPCs_id) <- c("id", paste0("diet", colnames(dietPCs$x)))

# Save diet PCA results as .rda (all outputs) & csv (of factor loadings) -------
saveRDS(dietPCs, file = paste0("../data/processed/ukb_dietpcs_24hr_adjkcal.rda"))
left_join(diet_all, vars_for_pca %>% mutate(id=as.integer(id)), by="id") %>%
  write.csv(file = paste0("../data/processed/ukb_diet24hr_dietpcs_adjkcal.csv"))



###############################################
## WaterFall plot of diet PC Factor Loadings ##
###############################################

# Descriptive diet labels
diet_labels_descriptive <- c(
  pc_vegetables = "Vegetables",
  pc_potatoes = "Potatoes",
  pc_fruit = "Fruit",
  pc_fruit_dried = "Dried fruit",
  pc_legumes = "Legumes",
  pc_whgrain = "Whole grains",
  pc_refgrain = "Refined grains",
  pc_dairy_fullfat = "Full-fat dairy",
  pc_dairy_lowfat = "Low-fat dairy",
  pc_nondairy_milk = "Non-dairy milks",
  pc_butter_margarine = "Butter & margarines",
  pc_eggs = "Eggs",
  pc_redmeat = "Red & processed meat",
  pc_poultry = "Poultry",
  pc_fish = "Fish & seafood",
  pc_plantbased_meat = "Plant-based meat alternatives",
  pc_nuts_seeds = "Nuts & seeds",
  pc_oliveoil = "Oils",
  pc_condiments = "Condiments",
  pc_added_sugar = "Added sugar",
  pc_sweets_desserts = "Sweets & desserts",
  pc_mixed_dishes = "Mixed dishes",
  pc_coffee = "Coffee",
  pc_tea = "Tea",
  pc_ssbs_sodas = "SSBs",
  pc_fruitjuice = "Fruit juices",
  pc_water = "Water",
  pc_beer = "Beer",
  pc_spirits = "Spirits",
  pc_wine = "Wine"
  )

#diet_labels_descriptive <- c(
#  fg_bevs_alcoholic = "Alcoholic beverages",
#  fg_cereals = "Cereals",
#  fg_dairy_products = "Dairy products",
#  fg_eggs = "Eggs",
#  fg_fats_spreads = "Fats and spreads",
#  fg_fish_dishes = "Fish dishes",
#  fg_fruits = "Fruits",
#  fg_meat_products = "Meat products",
#  fg_meat_substitutes = "Meat substitutes",
#  fg_mixed_dishes = "Mixed dishes",
#  fg_bevs_nonalch = "Non-alcoholic beverages",
#  fg_nuts_seeds = "Nuts & seeds",
#  fg_condiments = "Condiments",
#  fg_sweets_snacks = "Sweets & snacks",
#  fg_vegetables_potatoes = "Vegetables & potatoes"
#)


palette_waterfall = c("#888363","#C5C1A5", "#96A0B3", "#435269")


# =========================
## Plot dietPC waterfall
# =========================

## Run on local R

dietpcs <- readRDS("../data/processed/ukb_dietpcs_24hr.rda")

# prepare data for visualization
dietPC.loadings <- dietPCs$rotation %>% 
  as.data.frame() %>%
  mutate(Diet=diet_labels_descriptive[gsub("_adjkcal", "", rownames(.))], .before=PC1) %>%
  arrange(PC1)
dietPC.mat <- as.matrix(dietPC.loadings[,-1])
rownames(dietPC.mat) <- dietPC.loadings$Diet

plot_dietPC_waterfall.fun <- function(nPCs=10) {
  dietPC.loadings %>%
    select(Diet, c(paste0("PC", 1:nPCs))) %>%
    mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
    pivot_longer(-Diet) %>% 
    mutate(name=factor(name, levels=c(paste0("PC", 1:nPCs)), labels=c(paste0("Diet PC", 1:nPCs)) )) %>%
    mutate(direction = factor(ifelse(value>0.2, "Positive/Major", ifelse(value<0.2 & value>0, "Positive/Minor", 
                                                                         ifelse(value<0 & value>-0.2, "Negative/Minor", "Negative/Major"))),
                              levels=c("Positive/Major", "Positive/Minor", "Negative/Minor", "Negative/Major"))) %>%
    ggplot(aes(x=value, y=Diet, fill=direction)) + 
    facet_wrap(~name, ncol=nPCs/2) + 
    geom_col() + ylab("") + xlab("Rotated Factor Loadings") +
    geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values=palette_waterfall, name = "Factor Loadings") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          axis.text = element_text(color="black"),
          legend.position = "top",
          strip.text = element_text(face="bold", size=16),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          axis.title.x = element_text(size=14),
          legend.key.size = unit(1,"cm"),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16, face="bold"))
}

plot_dietPC_waterfall.fun(10)
