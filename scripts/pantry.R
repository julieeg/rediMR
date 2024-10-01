#Pantry file


# load packages
library(tidyverse) ; library(table1)


################################################################################
## Build & store standard parameter inputs & variable lists
################################################################################

# =========================================
##  Covariate sets for base GWAS
# =========================================
## Base covariates ===========
covarSetsBase <- list(
  gwas = list(
    Label="Base covariates",
    Covars=c("age","sex", paste0("gPC", 1:10)),
    Formatted=paste0(c("age","sex", paste0("gPC", 1:10)), collapse = "+")),
  agesex = list(
    Label="Base covariates",
    Covars=c("age","sex"),
    Formatted=paste0(c("age","sex"), collapse = "+")),
  dietpcs = list(
    Label="All Diet PCs",
    Covars=c(paste0("dietPC", 1:23)),
    Formatted=paste0("dietPC", 1:23, collapse="+"))
)

# ==============================================
##  Covariate sets for confounder adjustment
# ==============================================

covarSets <- list()

## Top 1-23 diet PCs (no confounders) ===========
for (i in 1:23) {
  covarSets[[i]] <- list(
    Label=paste0("Top", i, "DietPCs"),
    Covars=c(sapply(1:i, function(j) paste0("dietPC", j))),
    Names=c(sapply(1:i, function(j) paste0("Diet Pattern PC", j)))
  ) ; names(covarSets[[i]]$Names) <- c(covarSets[[i]]$Covars)
} ; names(covarSets) <- c("dietpc1", paste0("dietpctop", 2:23))

covarSets$alldietpcs = list(
  Label = "DietPCs",
  Covars = paste0(covarSets$dietpctop23$Covars, collapse = "+"),
  Names = "DietPCs"
)


## Confounders ==============

# Pilot ("stnd") covariates -------------------------
confounder_Label <- c("Smoking"="smoke", "Alcohol"="alch", "Physical Activity"="pa", 
                      "Income"="inc", "Education"="educ", 
                      "BMI"="bmi", "Waist2Hip"="w2h")

covarSets$basic = list(
  Label = "Basic covariates",
  Covars = c(
    "smoke_level.lab", "alch_freq.lab", "physact_level.lab", 
    "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10)),
  Names = c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist-to-hip",
    dietPC1="Diet Pattern PC1", dietPC2="Diet Pattern PC2", dietPC3="Diet Pattern PC3",
    dietPC4="Diet Pattern PC4", dietPC5="Diet Pattern PC5", 
    dietPC6="Diet Pattern PC6", dietPC7="Diet Pattern PC7", dietPC8="Diet Pattern PC8", 
    dietPC9="Diet Pattern PC9", dietPC10="Diet Pattern PC10") 
)

covarSets$basic_noalch = list(
  Label = "Basic covariates",
  Covars = c(
    "smoke_level.lab", "physact_level.lab", 
    "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10)),
  Names = c(
    smoke_level.lab="Smoking",
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist-to-hip",
    dietPC1="Diet Pattern PC1", dietPC2="Diet Pattern PC2", dietPC3="Diet Pattern PC3",
    dietPC4="Diet Pattern PC4", dietPC5="Diet Pattern PC5", 
    dietPC6="Diet Pattern PC6", dietPC7="Diet Pattern PC7", dietPC8="Diet Pattern PC8", 
    dietPC9="Diet Pattern PC9", dietPC10="Diet Pattern PC10") 
)


# Confounders 1 ---------------------------

covarSets$confounders1 = list(
  Label = "All Confounders",
  Covars = c(
    "smoke_level.lab", "alch_freq.lab", "physact_level.lab", 
    "income_level.lab", "educ_level.lab", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
)

covarSets$confounders1_num = list(
  Label = "All Confounders (numeric)",
  Covars = c(
    "smoke_level.num", "alch_freq.num", "physact_level", 
    "income_level.num", "educ_level.num", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
)

# Confounders 1 (no alcohol) ---------------------------

covarSets$confounders1_noalch = list(
  Label = "All Confounders",
  Covars = c(
    "smoke_level.lab", "physact_level.lab", 
    "income_level.lab", "educ_level.lab", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking", 
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
)

covarSets$confounders1_nolach = list(
  Label = "All Confounders (numeric)",
  Covars = c(
    "smoke_level.num", "physact_level", 
    "income_level.num", "educ_level.num", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking",
    physact_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
)


## Top 1-23 diet PCs and EACH confounder ===================
dietEachConf = list(
  Label="Diet PCs + Each Confounder",
  Covars=c(sapply(1:length(confounder_Label), function(i) {
    paste0(paste0(covarSets$dietpctop23$Covars, collapse="+"), "+", covarSets$confounders_num$Covars[i])})),
  Names=c(sapply(1:length(confounder_Label), function(i) paste0("Diet PCs + ", names(confounder_Label)[i])))
) ; names(dietEachConf$Names) <- paste0("dietpcsAnd", confounder_Label)

add_covarset2 <- list()
for(i in 1:length(confounder_Label)) {
  add_covarset2[[i]] <- list(
    Label = dietEachConf$Names[[i]],
    Covars = dietEachConf$Covars[i],
    Names = dietEachConf$Names[i])
} ; names(add_covarset2) <- names(dietEachConf$Names)


## Top 1-23 diet PCs with sequentially ADDED confounders =========================
dietAddConf = list(
  Label = "Diet PCs Adding Confounders", 
  Covars = c(sapply(1:length(confounder_Label), function(i) {
    paste0(covarSets$alldietpcs$Covars, "+", paste0(covarSets$confounders_num$Covars[1:i], collapse = "+"))})),
  Names = c(sapply(1:length(confounder_Label), function(i) paste0("Diet PCs + ", i, " Confounders (+", names(confounder_Label)[i], ")") ))
) ; names(dietAddConf$Names) <-  paste0("dietpcsAdd", confounder_Label)

add_covarset3 <- list()
for(i in 1:length(confounder_Label)) {
  add_covarset3[[i]] <- list(
    Labels = dietAddConf$Names[i],
    Covars = dietAddConf$Covars[i],
    Names = dietAddConf$Names[i]) 
} ; names(add_covarset3) <- names(dietAddConf$Names)

## append
covarSets <- c(covarSets, add_covarset2, add_covarset3) 
names(covarSets)


###############################
## Groups of covariate sets
###############################

covarSetGroups <- list(
  
  dietGroup1 = list(
    Sets = c("dietpc1", "dietpctop5", "dietpctop10", "dietpctop15", "dietpctop20", "dietpcall"),
    Names = c(dietpc1="Diet PC1", dietpctop5="Top 5 Diet PCs", dietpctop10="Top 10 Diet PCs", 
                       dietpctop15 = "Top 15 Diet PCs", dietpctop20="Top 20 Diet PCs"),
    Labels = c("Diet PC1", "Top 5 Diet PCs", "Top 10 Diet PCs", "Top 15 Diet PCs", "Top 20 Diet PCs")),
  
  dietGroup2 = list(
    Sets = c(names(covarSets)[1:23]),
    Names = c(sapply(c(1:23), function(i) covarSets[[i]]$Names)),
    Labels =  c(sapply(c(1:23), function(i) covarSets[[i]]$Label))
    ),
  
  dietEachConf = list(
    Sets = c(names(dietEachConf$Names)),
    Names = c(dietEachConf$Names),
    Labels = as.vector(dietEachConf$Names)),
  
  dietAddConf = list(
    Sets = c(names(dietAddConf$Names)),
    Names = c(dietAddConf$Names),
    Labels = as.vector(dietAddConf$Names))
)
  
# ==============================
## Descriptive variable labels
# ==============================

diet_traits <- c(
  oilyfish = "oily fish intake (sev/wk)",
  bread_intake = "bread type (white vs. whole/brown)",
  alch_glasspermonth = "alcohol intake (glass/month)"
)

outcome_traits <- c(
  tg="Triglyceride",
  cvd="CVD",
  ldl="LDL",
  cir="Cirrhosis",
  alt="ALT"
)

confounders1 <- c(
  smoke_level.lab = "Smoking",
  alch_freq.lab = "Alcohol",
  physact_level.lab = "Physical Activity",
  income_level.lab = "Income",
  educ_level.lab = "Education",
  bmi = "BMI", 
  waist2hip = "Waist-to-hip"
)

confounders1_num <- c(
  smoke_level.num = "Smoking",
  alch_freq.num = "Alcohol",
  physact_level.num = "Physical Activity",
  income_level.num = "Income",
  educ_level.num = "Education",
  bmi="BMI",
  waist2hip = "Waist-to-hip"
)


# ======================================
## Compile MR outcomes from ieuGWAS
# ======================================

outcomeVars <- list(
  GCST90239664=list(
    Label="log(TG)",
    description="log(Triglyceride)",
    gwas="Graham S, 2021",
    var="tg_log"),
  GCST90132314=list(
    Label="CAD",
    description="Coronary Artery Disease",
    gwas="Aragam KG, 2022",
    var="cvd"),
  GCST90239658=list(
    Label="LDL",
    description="LDL Cholesterol",
    gwas="Graham S, 2021",
    var="ldl"),
  GCST90319877=list(
    Label="Cirrhosis",
    description="Cirrhosis of Liver",
    gwas="Ghouse J, 2024",
    var="cir"),
  GCST90013405=list(
    Label="ALT",
    description="Alanine Aminotransferase",
    gwas="Pazoki R, 2021",
    var="alt")
)

# ======================================
## Diet-trait pairs
# ======================================

diet_trait_pairs=c("alch_alt","alch_cir", "breadtype_cvd", "breadtype_ldl", "oilyfish_cvd","oilyfish_tg")
pairs.l <- list(
  oilyfish_tg = list(
    exposure="oilyfish_QT",
    outcome="tg",
    pair.label="Oily fish intake on TG",
    exposure.label="Oily fish intake",
    outcome.label="Triglyceride"),
  oilyfish_cvd = list(
    exposure="oilyfish_QT",
    outcome="cvd",
    pair.label="Oily fish intake on CVD",
    exposure.label="Oily fish intake",
    outcome.label="Coronary Artery Disease"),
  bread_ldl = list(
    exposure="bread_type_BIN",
    outcome="ldl",
    pair.label="Bread type (choose whole grain) on LDL",
    exposure.label="Bread type",
    outcome.label="LDL Cholesterol"),
  bread_cvd = list(
    exposure="bread_type_BIN",
    outcome="cvd",
    pair.label="Bread type (choose whole grain) on CVD",
    exposure.label="Bread type",
    outcome.label="Coronary Artery Disease"),
  alch_cir = list(
    exposure="alch_glasspermonth_QT",
    outcome="cir",
    pair.label="Alcohol intake (glass/month) on Cirrhosis",
    exposure.label="Alcohol (glasses/m)",
    outcome.label="Cirrhosis"),
  alch_alt = list(
    exposure="alch_glasspermonth_QT",
    outcome="alt",
    pair.label="Alcohol intake (glass/month) on ALT",
    exposure.label="Alcohol (glasses/m)",
    outcome.label="Alanine Aminotransferase")
)

sets <- c("all"="All", "refined_lt20"="Refined <20%", 
          "refined_lt10"="Refined <10%", "refined_lt05"="Refined <5%")

################################################################################
## Basic Functions
################################################################################

#####################################################
## ~~~~  Data display for descriptive tables  ~~~~ ##
#####################################################

# ======================================
## Print continuous vars as mean +- SD 
# ======================================
mean_sd<-function(x, d=2) {
  sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
          round(sd(x, na.rm=T), digits = d))
}


# ==================================
## Print categorical vars as n (%)
# ==================================
n_pct <- function(x, level=F) {
  if(level==F) {
    sapply(as.list(names(table(x))), function(lvl) {
      paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100,1), "%)") }) } 
  else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100,1), "%)")}
}


# =====================================================
## Print continuous vars as median [25th, 75th %tile]
# =====================================================
median_25to75<-function(x, d=2) {
  qs<-round(quantile(x, breaks=seq(0,1,0.25), na.rm=T), d)
  sprintf("%s [%s, %s]", round(median(x, na.rm=T), digits = d), 
          qs[[2]], qs[[4]])
}



################################################
##  ~~~~ Data wrangling & transformation ~~~~ ##
################################################

# ========================
## Remove outliers by SD
# ========================
remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}


# =================================================
## Median impute for negative or missing values
# =================================================
median_imp.fun <- function(x) {
  x.new <- ifelse(x == -1 | x == -3 | x == -9 | is.na(x) == T, median(x, na.rm=T), x)
  return(x.new)
}


# ==========================
## Winsorize data by SD
# ==========================
winsorize <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x<bounds[1], bounds[1], ifelse(x>bounds[2], bounds[2], x))
  x
}


# ========================
## Recode as zscore
# ========================
zscore.fun <- function(x) {
  z<-((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  return(z)
}


# ===========================================
## Add descriptive labels to UKB variables
# ===========================================
descr_label.fun <- function(data, base_var, labs_vals) {
  
  base <- data %>% select(base_var) 
  temp <- rep(NA, length(base))
  
  for(i in 1:length(labs_vals)) {
    temp[base == labs_vals[[i]] ] <- names(labs_vals)[i]
  } ; return(temp)
}



################################################################
## Descriptive "Table 1" functions
################################################################

# =========================
## Trait-confounder associations 
# ==================

corPhenoConf.fun <- function(pheno, conf_num, adjCovar, data=dat) {
  
  if(length(adjCovar) > 1) {adjCovar =paste0(adjCovar, collapse="+")}
  
  base=summary(lm(formula(paste0(pheno, "~", conf_num)), data=dat))
  adj <- lm(formula(paste0(pheno, "~", snp_EA, "+", gwasCovarsBase, "+", adjCovar)), data)
  
  baseP=summary(baseM)$coef[2,4]
  adjP=summary(adjM)$coef[2,4]
  
  baseSE=summary(baseM)$coef[2,2]
  adjSE=summary(adjM)$coef[2,2]
  
  baseCI <- confint(baseM)[2,]
  adjCI <- confint(adjM)[2,]
  pctB <-round(((adjM$coef[2]-baseM$coef[2])/baseM$coef[2])*100, 2)
  
  
  if(!is.na(covarName)) {Covar = covarName} else {Covar = adjCovar}
  lmOut <- cbind.data.frame(
    snp=snp_EA, Covar=covarName, B_base=baseM$coef[2], B_adj=adjM$coef[2], 
    SE_base=baseSE, SE_adj=adjSE, lowCI_base=baseCI[1], upCI_base=baseCI[2], 
    lowCI_adj=adjCI[1], upCI_adj=adjCI[2], P_base=baseP, P_adj=adjP, 
    B_pctChange=pctB) ; rownames(lmOut) <- paste0(Covar)
  
  return(lmOut)
}

# ================================================
## Reformatting Bdat data to long for plotting
# ================================================

pivot_Bdat_to_long <- function(Bchange.df) {
  x <- Bchange.df %>%
    pivot_longer(c(ends_with("base"), ends_with("adj")), names_to=c("stat", "mod"), names_sep="_") %>%
    pivot_wider(names_from="stat", values_from = value) %>%
    mutate(Covar = ifelse(Covar == "Smoking" & mod == "base", "Unadjusted", Covar)) %>%
    filter(Covar == "Unadjusted" | mod == "adj") %>%
    mutate(B_pctChange = ifelse(Covar == "Unadjusted", 0, B_pctChange)) %>%
    mutate(Covar = ifelse(Covar == "Top23DietPCs", "DietPCs", Covar))
  return(x)
}


######################################
## Primary rediMR functions
######################################

# =================================================
## Derive diet PCs without exposure diet trait
# =================================================

derive_dietPCs <- function(exposure, exclude, data=dat) {
  
  vars_for_pca <- data %>% select(
    id,
    cooked_veg_QT=cooked_veg, raw_veg_QT=raw_veg,
    fresh_fruit_QT=fresh_fruit, dried_fruit_QT=dried_fruit, 
    oily_fish_QT=oily_fish, nonoily_fish_QT=nonoily_fish,
    procmeat_QT=procmeat, poultry_QT=poultry, cheese_QT=cheese,
    beef_QT=beef, lamb_QT=lamb, pork_QT=pork, 
    bread_type_white_vs_brown_or_whole_BIN=bread_type_white_vs_brown_or_whole, 
    bread_intake_QT=bread_intake,
    milk_type_full_vs_low_or_nonfat_BIN=milk_type_full_vs_low_or_nonfat,
    cereal_type_sugar_vs_any_bran_BIN=cereal_type_sugar_vs_any_bran, 
    cereal_intake_QT=cereal_intake,
    spread_type_butter_vs_any_other_BIN=spread_type_butter_vs_any_other,
    coffee_type_decaf_vs_regular_BIN=coffee_type_decaf_vs_regular, coffee_QT=coffee,
    tea_QT=tea, water_QT=water, 
    addsalt_always_often_vs_nrs_BIN=addsalt_always_often_vs_nrs,
    hotdrink_temp_hot_or_vhot_vs_warm_BIN=hotdrink_temp_hot_or_vhot_vs_warm) %>%
    
    # Replace missing data with median values
    mutate_at(vars(-id), function(x) ifelse(is.na(x), median(x, na.rm=T), x)) %>%
    
    # Winzorise
    mutate(across(where(is.numeric), ~winsorize(., SDs=5))) %>%
    filter(complete.cases(.)==T)
  
  if (exclude == "none") {
    cat("Building diet patterns via pca with all diet traits ...")
    vars_for_pca <- vars_for_pca } else {
      cat(paste0("Building diet patterns via pca, excluding ", exclude))
      vars_for_pca <- vars_for_pca %>% select(-(starts_with(exclude)))
    }
  
  ## Run PCA 
  diet_pcs <- prcomp(select(vars_for_pca, -id), scale.=T)  # Run PCA
  
  ## Extract scores
  diet_pcs.scores <- cbind.data.frame(vars_for_pca$id, diet_pcs$x)
  names(diet_pcs.scores) <- c("id", paste0("dietPC", 1:ncol(diet_pcs$x)))
  
  return(list(pcs=diet_pcs, scores=diet_pcs.scores))
  
}
  

# =================================================
## Calculate %B change with covariate adjustment
# =================================================

pctBchange.fun <- function(pheno, snp, adjCovar, baseCovars=covarSetsBase$gwas$Covars, covarName=NA, data=dat) {
  
  if(length(baseCovars) > 1) {baseCovarsFormat <- paste0(baseCovars, collapse="+")} else {
    baseCovarsFormat <- baseCovars} 
  if(length(adjCovar) > 1) {adjCovar <- paste0(adjCovar, collapse="+")}
  
  snp_EA <- names(data %>% select(starts_with(snp)))
  baseM <- lm(formula(paste0(pheno, "~", snp_EA, "+", baseCovarsFormat)), data)
  adjM <- lm(formula(paste0(pheno, "~", snp_EA, "+", baseCovarsFormat, "+", adjCovar)), data)
  
  baseP=summary(baseM)$coef[2,4]
  adjP=summary(adjM)$coef[2,4]
  
  baseSE=summary(baseM)$coef[2,2]
  adjSE=summary(adjM)$coef[2,2]
  
  baseCI <- confint(baseM)[2,]
  adjCI <- confint(adjM)[2,]
  pctB <-round(((adjM$coef[2]-baseM$coef[2])/baseM$coef[2])*100, 2)
  
  
  if(!is.na(covarName)) {Covar = covarName} else {Covar = adjCovar}
  lmOut <- cbind.data.frame(
    snp=snp_EA, Covar=covarName, B_base=baseM$coef[2], B_adj=adjM$coef[2], 
    SE_base=baseSE, SE_adj=adjSE, lowCI_base=baseCI[1], upCI_base=baseCI[2], 
    lowCI_adj=adjCI[1], upCI_adj=adjCI[2], P_base=baseP, P_adj=adjP, 
    B_pctChange=pctB) ; rownames(lmOut) <- paste0(Covar)
  
  return(lmOut)
}



#########################################################################
##  ~~~~  Basic descriptive plots of continuous/categorical vars  ~~~~ ##
#########################################################################

# ============================================
## pre-built color palettes 
# ============================================

library("paletteer") ; library("RColorBrewer")
# Developmental version
#devtools::install_github("awhstin/awtools")

palettes <- list(NatComms= paletteer_d("ggsci::nrc_npg", n=10),
                 classicgreen=paletteer::paletteer_c("ggthemes::Classic Green", 100)[seq(1,100,4)],
                 classicorange=paletteer::paletteer_c("ggthemes::Classic Orange", 100)[seq(1,100,4)],
                 greens5 = paletteer_dynamic("cartography::green.pal", 5),
                 greens = rev(paletteer_dynamic("cartography::green.pal", 10)),
                 pugr = c("#9C50CE", "#9C50CE95","#30900195", "#309001"),
                 oranges = rev(paletteer_dynamic("cartography::orange.pal", 10)),
                 blues = rev(paletteer_dynamic("cartography::blue.pal", 10)),
                 purples=rev(paletteer_dynamic("cartography::purple.pal", 10))[3:10],
                 Conf=c(brewer.pal(9,"Oranges")[4:6], brewer.pal(9, "Blues")[5:6],  
                            brewer.pal(9, "Purples")[c(8:2)]),
                            #brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6]),
                 ConfNames=c( 
                   "Unadjusted"="#000000",
                   "Smoking"=brewer.pal(9,"Oranges")[4],
                   "Alcohol"=brewer.pal(9,"Oranges")[5],
                   "PA Level"=brewer.pal(9,"Oranges")[6],
                   "Income"= brewer.pal(9, "Blues")[5],
                   "Education"= brewer.pal(9, "Blues")[6],
                   "BMI"=brewer.pal(9, "Purples")[8],
                   "Waist-hip"=brewer.pal(9, "Purples")[7]),
                 #brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6]),
                 Confnum=c(
                   "smoke_level.num" = brewer.pal(9, "Oranges")[4], 
                   "alch_freq.num" = brewer.pal(9, "Oranges")[5], 
                   "physact_level" = brewer.pal(9, "Oranges")[6],
                   "income_level.num" = brewer.pal(9, "Blues")[5],
                   "educ_level.num" = brewer.pal(9, "Blues")[6],
                   "bmi" = brewer.pal(9, "Purples")[8],
                   "waist2hip" = brewer.pal(9, "Purples")[7]
                 ),
                 CovarsDiet = c("Diet Pattern PC1" = brewer.pal(9, "Purples")[6],
                                "Diet Pattern PC2" = brewer.pal(9, "Purples")[5],
                                "Diet Pattern PC3" = brewer.pal(9, "Purples")[4],
                                "Diet Pattern PC4" = brewer.pal(9, "Purples")[3],
                                "Diet Pattern PC5" = brewer.pal(9, "Purples")[2],
                                "Diet Pattern PC6" = brewer.pal(9, "PiYG")[1],
                                "Diet Pattern PC7" = brewer.pal(9, "PiYG")[2],
                                "Diet Pattern PC8" = brewer.pal(9, "PiYG")[2],
                                "Diet Pattern PC9" = brewer.pal(9, "PuRd")[5],
                                "Diet Pattern PC10" = brewer.pal(9, "PuRd")[6]),
                 hm_bwr = colorRampPalette(c("Blue", "White", "Red"))(1024)
)



# ========================
## Build ggplot themes
# ========================

ggtheme <- theme_bw() + theme(panel.grid.minor.y = element_blank(), 
                 panel.grid.minor.x = element_blank(), 
                 axis.text.x = element_text(color="black", size=8), #angle=35, hjust=1, 
                 axis.text.y = element_text(color="black", size=8),
                 axis.title = element_text(color="black", size=8),
                 legend.position = "right", 
                 legend.box.background = element_rect(color = "black"),
                 legend.key.size = unit(0.25, 'line'),
                 legend.margin = margin(0.2,0.2,0.2,0.2, unit="pt"),
                 legend.text = element_text(size=8), 
                 legend.title = element_text(face="bold", size=8),
                 plot.title=element_text(size=8),
                 strip.text = element_text(face="bold", size=8)
)



ggtheme2 <- theme_bw() + theme(panel.grid.minor.y = element_blank(), 
                              panel.grid.minor.x = element_blank(), 
                             # panel.grid.major.y = element_blank(), 
                              panel.grid.major.x = element_blank(),
                              axis.text.x = element_text(color="black", size=9), #angle=35, hjust=1, 
                              axis.text.y = element_text(color="black", size=9),
                              axis.title = element_text(color="black", size=8, face="bold"),
                              legend.position = "right", 
                              legend.box.background = element_rect(color = "black"),
                              legend.key.size = unit(0.25, 'line'),
                              legend.margin = margin(0.2,0.2,0.2,0.2, unit="pt"),
                              legend.text = element_text(size=8), 
                              legend.title = element_text(face="bold", size=9),
                              plot.title=element_text(size=9),
                              strip.text = element_text(face="bold", size=9)
)


ggthemeMR <- theme_bw() + 
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(color="black", size=9), #angle=35, hjust=1, 
        axis.text.y = element_text(color="black", size=9),
        axis.title = element_text(color="black", size=10, face="bold"),
        legend.position = "none", 
        #legend.box.background = element_rect(color = "black"),
        #legend.key.size = unit(0.25, 'line'),
        #legend.margin = margin(0.2,0.2,0.2,0.2, unit="pt"),
        #legend.text = element_text(size=8), 
        #legend.title = element_text(face="bold", size=9),
        #plot.title=element_text(size=9),
        strip.text = element_text(face="bold", size=10)
  )
  
  
ggthemeMRksj <- theme_bw() +
  ggplot2::theme(
    axis.title.y=ggplot2::element_text(size=10), #14
    axis.text.y=ggplot2::element_text(size=10), 
    axis.ticks.y=ggplot2::element_line(linewidth=0),
    axis.title = ggplot2::element_text(color="black", size=10, face="bold"),
    axis.text.x=ggplot2::element_text(size=10),
    axis.ticks.x=ggplot2::element_line(linewidth=0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2))

# =========================================
## Plot continuous x categorical
# =========================================

# Functions 
plot_contXcat <- function(cat_var, cont_var=pheno, data=dat, ...) {
  data %>% select(cat=all_of(cat_var), cont=all_of(cont_var)) %>% 
    filter(!is.na(cont)) %>%
    group_by(cat) %>% summarise(m=mean(cont, na.rm=T), sd=sd(cont, na.rm=T)) %>%
    ggplot(aes(x=cat, y=m, ymin=m-sd, ymax=m+sd, fill=cat)) + ggthemeF + 
    theme(axis.text.x=element_text(angle=35, hjust=1, color="black")) +
    geom_bar(stat="identity") + geom_errorbar(width=0.15) +
    scale_fill_manual(...) +
    labs(title=paste("Bar plot of", cont_var, "by", cat_var), x=" ", y="mean (SD)")
}


plot_contXcont <- function(cont_var_x, cont_var_y=pheno, data=dat, ...) {
  data %>% select(cont_x=all_of(cont_var_x), cont_y=all_of(cont_var_y)) %>% filter(complete.cases(.)) %>%
    ggplot(aes(x=cont_x, y=cont_y)) + ggthemeF +
    geom_point(size=1.25, position=position_jitter(0.15), alpha=0.75, ...) + 
    geom_smooth(method="lm", color = "black") + 
    labs(title=paste("Scatter plot of", cont_var_y, "by", cont_var_x), y=cont_var_y, x=cont_var_x)
}


plot_catxcat_pct <- function(cat_var_group, cat_var_x="phenoQ", data=dat, ...) {
  prop.table(table(data %>% select(cat_x=all_of(cat_var_x), cat_group=all_of(cat_var_group))),1) %>% 
    as.data.frame() %>%
    ggplot(aes(x=cat_x, y=Freq*100, group=cat_group, fill=cat_group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title=paste(cat_var_x, "by", cat_var_group), y="Percent", x=" ") + 
    scale_fill_manual(...) +
    theme(plot.title = element_text(face="bold", size=10),
          axis.text.x = element_text(angle=25, vjust=0.95)) +
    ggthemeF 
}


# =============
## plot dietPC PVE or eigenvalues 
# ============

plot_dietPC_pve.fun <- function(dPC_sdev, title) {
  as.data.frame(cbind(nPC=1:length(dPC_sdev),
                      pve=((dPC_sdev^2)/sum(dPC_sdev^2))*100)) %>%
  ggplot(aes(x=nPC, y=pve)) +
    geom_point(size=2, shape=19, alpha=0.75) + 
    geom_line(linewidth=0.35) + 
    scale_x_continuous(breaks=seq(1,length(dPC_sdev),1)) + 
    xlab("diet PC #") + ylab("Percent Variance Explained") + 
    scale_y_continuous(limits=c(2,11), breaks=seq(2,12,2)) +
    ggtitle(title) + ggtheme2 
}

plot_dietPC_egn.fun <- function(dPC_sdev, title) {
  as.data.frame(cbind(nPC=1:length(dPC_sdev),
                      egn=(dPC_sdev^2))) %>%
    ggplot(aes(x=nPC, y=egn)) +
    geom_point(size=2, shape=19, alpha=0.75) + 
    geom_line(linewidth=0.35) + 
    scale_x_continuous(breaks=seq(1,length(dPC_sdev),1)) + 
    xlab("diet PC #") + ylab("Eigenvalues") +
    scale_y_continuous(limits=c(0.425,2.5), breaks=seq(0.5,2.5,0.5)) +
    geom_hline(yintercept = 1, linetype="dashed", linewidth=0.35) + 
    ggtitle(title) + ggtheme2 
}



#EOF


