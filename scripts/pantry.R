#Pantry file


# load packages
library(tidyverse) ; library(table1)


################################################################################
# ==============================================================================
## Build & store standard parameter inputs & variable lists
# ==============================================================================
################################################################################

# =========================================
##  Covariate sets
# =========================================

## Base covariates ===========
covarSetsBase <- list(
  
  gwas = list(
    Label="Base covariates",
    Covars=c("age","sex", paste0("gPC", 1:10)),
    Formatted=paste0(c("age","sex", paste0("gPC", 1:10)), collapse = "+")
    ),
  
  agesex = list(
    Label="Base covariates",
    Covars=c("age","sex"),
    Formatted=paste0(c("age","sex"), collapse = "+")),
  
  dietpcs = list(
    Label="All Diet PCs",
    Covars=c(paste0("dietPC", 1:23)),
    Formatted=paste0("dietPC", 1:23, collapse="+"))
  
)


## Top 1-23 diet PCs (no confounders) ===========
covarSets <- list()
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

confounder_Label <- c("Smoking"="smoke", "Alcohol"="alch", "Physical Activity"="pa", 
                      "Income"="inc", "Education"="educ", 
                      "BMI"="bmi", "Waist2Hip"="w2h")

covarSets$confounders = list(
  Label = "All_Confounders",
  Covars = c(
    "smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab", 
    "income_level.lab", "educ_level.lab", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    pa_met_excess_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
  )

covarSets$confounders_num = list(
  Label = "All_Confounders (numeric)",
  Covars = c(
    "smoke_level.num", "alch_freq.num", "pa_met_excess_level.num", 
    "income_level.num", "educ_level.num", "bmi", "waist2hip"),
  Names = c(
    smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
    pa_met_excess_level.lab="Physical Activity", income_level.lab = "Income", 
    educ_level.lab="Education", bmi="BMI", waist2hip="Waist2Hip")
)

## genetic PCs 
add_covarset <- list()
for(i in 1:20) {
  add_covarset[[i]] <- list(
    Label = paste0("Top", i, "geneticPCs"),
    Covars = c(sapply(1:i, function(j) paste0("gPC", j))),
    Names = c(sapply(1:i, function(j) paste0("Genetic PC", j)))
    ) ; names(add_covarset[[i]]$Names) <- c(add_covarset[[i]]$Covars)
} ; names(add_covarset) <- c(paste0("gPCtop", 1:20))

## Top 1-23 diet PCs and EACH confounder ===================
dietEachConf = list(
    Label="Diet PCs + Each Confounder",
    Covars=c(sapply(1:length(confounder_Label), function(i) {paste0(paste0(covarSets$dietpctop23$Covars, collapse="+"), "+", covarSets$confounders_num$Covars[i])})),
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
  Covars = c(sapply(1:length(confounder_Label), function(i) {paste0(covarSets$alldietpcs$Covars, "+", paste0(covarSets$confounders_num$Covars[1:i], collapse = "+"))})),
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
covarSets <- c(covarSets, add_covarset, add_covarset2, add_covarset3) 

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
  
  dietPCs = list(
    Sets = c(names(covarSets)[1:23]),
    Names = c(sapply(c(1:23), function(i) covarSets[[i]]$Names)),
    Labels = c(sapply(c(1:23), function(i) covarSets[[i]]$Label))),
  
  dietEachConf = list(
    Sets = c(names(dietEachConf$Names)),
    Names = c(dietEachConf$Names),
    Labels = as.vector(dietEachConf$Names)),
  
  dietAddConf = list(
    Sets = c(names(dietAddConf$Names)),
    Names = c(dietAddConf$Names),
    Labels = as.vector(dietAddConf$Names)),
  
  geneticPCs = list(
    Sets = c(names(covarSets)[27:46]),
    Names = c(sapply(c(27:46), function(i) covarSets[[i]]$Names)),
    Labels = c(sapply(c(27:46), function(i) covarSets[[i]]$Label))
  )
)
  

# =========================
## Diet traits 
# =========================

diet_labels <- c(
  raw_veg="Raw vegetable intake",
  cooked_veg="Cooked vegetable intake",
  fresh_fruit="Fresh fruit intake",
  dried_fruit="Dried fruit intake",
  oily_fish="Oily fish intake",
  nonoily_fish="Non-oily fish intake",
  procmeat="Processed meat intake",
  poultry="Poultry intake",
  cheese="Cheese intake",
  beef="Beef intake",
  lamb="Lamb intake",
  pork="Pork intake",
  bread_type_white_vs_brown_or_whole="Choose white > brown/whole bread",
  bread_intake="Bread intake",
  milk_type_full_vs_low_or_nonfat="Choose full>low/nonfat milk",
  cereal_type_sugar_vs_any_bran="Choose sugary>bran cereal",
  coffee_type_decaf_vs_regular="Choose decaf>regular coffee",
  cereal_intake="Cereal intake",
  spread_type_butter_vs_any_other="Choose butter>other spreads",
  coffee="Coffee intake",
  tea="Tea intake",
  water="Water intake",
  addsalt_always_often_vs_nrs="Add salt always/often",
  hotdrink_temp_hot_or_vhot_vs_warm="Prefer very/hot>warm drinks"
)

## vectors 
foods <- c("Raw vegetables"="raw_veg", "Fresh Fruit"="fresh_fruit")
foods.l <- as.list(foods)


# =========================
## Confounders
# =========================

confounders <- c(
  smoke_level.lab = "Smoking",
  alch_freq.lab = "Alcohol",
  pa_met_excess_level.lab = "Physical Activity",
  income_level.lab = "Income",
  educ_level.lab = "Education",
  bmi = "BMI", 
  waist2hip = "Waist-to-hip"
)

confounders_num <- c(
  smoke_level.num = "Smoking",
  alch_freq.num = "Alcohol",
  pa_met_excess_level.num = "Physical Activity",
  income_level.num = "Income",
  educ_level.num = "Education",
  bmi="BMI",
  waist2hip = "Waist-to-hip"
)


# ======================================
## Compile MR outcomes from ieuGWAS
# ======================================

outcomes.l <- list(
  colcan = list(
    Label="Colorectal Cancer", var="colorectal_cancer", id.outcome="finn-b-C3_COLORECTAL", type="OR",
    Name=c("colcan" = "Colorectal Cancer")),
  ischem = list(
    Label="Ischemic Stroke", var="ischemic_stroke", id.outcome="finn-b-I9_STR_EXH_EXNONE", type="OR",
    Name=c("ischem"="Ischemic Stroke")),
  depr = list(
    Label = "Depression", var="depression", id.outcome="finn-b-F5_DEPRESSIO", type="OR",
    "depr"= "Depression"),
  t2d = list(
    Label="T2D", var="t2d", id.outcome="ebi-a-GCST006867", type="OR", Name=c("t2d"= "T2D")),
  crp = list(
    Label="CRP", var="crp", id.outcome="prot-a-670", type="Beta", Name=c("crp"="CRP")),
  tg = list(
    Label="TG", var="tg", id.outcome="ebi-a-GCST002216", type="Beta", Name=c("tg"="TG")),
  dbp = list(
    Label="DBP", var="dbp", id.outcome="ebi-a-GCST90018952", type="Beta",Name=c("dbp"="DBP")),
  gluc = list(
    Label="Glucose", var="glucose", id.outcome="ieu-b-113", type="Beta", Name=c("gluc"="Glucose"))
)

outcome.Labs <- c(sapply(1:length(outcomes.l), function(i) outcomes.l[[i]]$label))
outcome_vars <- names(outcomes.l)



################################################################################
# ==============================================================================
## Basic Functions
# ==============================================================================
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

# ======================================
## Load print_summary_table function
# ======================================

source("../scripts/functions/print_summary_table_fun.R", echo=F)


# ==============================================
## Summary table (data format)
# ==============================================

summary_table.fun <- function(pheno, vars_to_summarize) {
  
  do.call(rbind.data.frame, lapply(1:length(vars_to_summarize), function(i) {
    
    covar=vars_to_summarize[i]
    tmp <- dat %>% select(Level=all_of(names(covar)[[1]]), y=all_of(pheno)) %>% filter(!is.na(y)) 
    
    if(is.factor(tmp$Level)) {
      out <- tmp %>% group_by(Level) %>% 
        summarise(est.corr = mean(y, na.rm=T),
                  est.err = sd(y, na.rm=T)) %>% 
        mutate(Estimate="m SD", Variable=covar,  .before=Level) %>%
        as.data.frame()
    } 
    else if(is.numeric(tmp$Level)) {
      out.cor <- cor.test(tmp$Level, tmp$y)
      out.ci <-  out.cor$conf.int
      out <- as.data.frame(matrix(
        c("cor se", names(covar), "-", out.cor$est[[1]], ((out.ci[[1]]-out.cor$est[[1]])/-1.96) ),
        1, 5, dim=list(NULL, c("Estimate", "Variable", "Level", "est.corr", 'est.err'))))
    }
  } ))
}

# ==============================================
## Summary table (manuscript/pretty format)
# ==============================================

## Print nice summary table
print_summary_table.fun <- function(pheno, vars_to_summarize) {
  
  do.call(rbind.data.frame, lapply(1:length(vars_to_summarize), function(i) {
    covar=vars_to_summarize[i]
    tmp <- dat %>% select(Level=all_of(names(covar)[[1]]), y=all_of(pheno)) %>% filter(!is.na(y)) 
    if(is.factor(tmp$Level)) {
      out <- tmp %>% group_by(Level) %>% 
        summarise(Summary = mean_sd(y)) %>% 
        mutate(Variable=covar, .before=Level) %>%
        mutate(Estimate="m SD", .before=Variable) %>%
        as.data.frame()
    } 
    else if(is.numeric(tmp$Level)) {
      out.cor <- round(cor.test(tmp$Level, tmp$y)$est, 2)
      out.ci <-  round(cor.test(tmp$Level, tmp$y)$conf.int, 2)
      out <- as.data.frame(matrix(c(
        "Cor (95%CI)", covar, "-",  sprintf("%s (%s, %s)", out.cor, out.ci[1], out.ci[2])),
        1, 4, dim=list(NULL, c("Estimate", "Variable", "Level", "Summary"))))
    }
  } ))
}


# ====================================
## Trait-confounder associations 
# ====================================

corPhenoConf.fun <- function(pheno, conf_num, adjCovar, data=dat) {
  
  if(length(adjCovar) > 1) {adjCovar =paste0(adjCovar, collapse="+")}
  
  base=summary(lm(formula(paste0(pheno, "~", conf_num)), data=dat))
  adj <- lm(formula(paste0(pheno, "~", snp_EA, "+", baseCovars, "+", adjCovar)), data)
  
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
## Calculate %B change with covariate adjustment
# =================================================

pctBchange.fun <- function(pheno, snp, adjCovar, baseCovars=covarSetsBase$gwas$Covars, covarName=NA, data=dat) {
  
  if(length(baseCovars) > 1) {baseCovarsFormat =paste0(baseCovars, collapse="+")}
  if(length(adjCovar) > 1) {adjCovar =paste0(adjCovar, collapse="+")}
  
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
                 hm_bwr = colorRampPalette(c("Blue", "White", "Red"))(1024),
                 Conf=c(brewer.pal(9,"Oranges")[4:6], brewer.pal(9, "Blues")[5:6],  
                            brewer.pal(9, "Purples")[c(8:2)]),
                            #brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6]),
                 ConfNames=c( 
                   "Unadjusted"="#000000",
                   "Smoking"=brewer.pal(9,"Oranges")[4],
                   "Alcohol"=brewer.pal(9,"Oranges")[5],
                   "Phys Act"=brewer.pal(9,"Oranges")[6],
                   "Income"= brewer.pal(9, "Blues")[5],
                   "Education"= brewer.pal(9, "Blues")[6],
                   "BMI"=brewer.pal(9, "Purples")[8],
                   "Waist2Hip"=brewer.pal(9, "Purples")[7]),
                 ConfNames_num=c(
                   "smoke_level.num" = brewer.pal(9, "Oranges")[4], 
                   "alch_freq.num" = brewer.pal(9, "Oranges")[5], 
                   "pa_met_excess_level.num" = brewer.pal(9, "Oranges")[6],
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
                                "Diet Pattern PC10" = brewer.pal(9, "PuRd")[6])
)



# ========================
## Build ggplot themes
# ========================

ggtheme1 <- theme_bw() + theme(panel.grid.minor.y = element_blank(), 
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

ggtheme2 <- theme_bw() + theme(
  panel.grid.minor.y = element_blank(), 
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


# =======================================
## plot dietPC PVE or eigenvalues 
# ======================================

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


# =======================================
## waterfall plot of Factor loadings
# ======================================
plot_dietPC.fun <- function(pheno, dietPCset, nPCs=10) {
  as.data.frame(dietPCset) %>%
    mutate(Diet=diet_labels[gsub("_QT", "", gsub("_BIN", "", rownames(.)))]) %>%
    select(Diet, c(paste0("PC", 1:nPCs))) %>%
    mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
    pivot_longer(-Diet) %>% 
    mutate(name=factor(name, levels=c(paste0("PC", 1:nPCs))),
           dir_col = factor(case_when(
             value <= -0.2~"Negative/Major", value > -0.2 & value <0 ~ "Negative/Minor",
             value > 0 & value <0.2 ~"Positive/Minor", value >= 0.2 ~ "Positive/Major"),
             levels=c("Negative/Major", "Negative/Minor", "Positive/Minor", "Positive/Major"))) %>%
    ggplot(aes(x=value, y=Diet, fill=dir_col)) + ggtheme1 + 
    facet_wrap(~name, nrow = 2) + geom_col() + ylab(" ") + 
    xlab(paste("Rotated Factor Loadings")) +
    geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values=palettes$pugr, name = "Factor Loading") +
    theme(legend.position = "bottom",
          axis.text.y=element_text(size=10)) +
    ggtitle(paste("PCA-deried dietary patterns without: ", diet_labels[[pheno]])) 
}



#EOF


