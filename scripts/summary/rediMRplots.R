# rediMR accompanying descriptive plots



################################
##  Set up & Load data files  ##
################################

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer", "ggpubr"),  
       library, character.only = TRUE)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
tag <- args[2]
pheno_tag <- paste0(pheno, "_", tag)


## create lists of diet data
rediMR_dir="../data/processed/rediMR"
out_dir="../output"


# load basic functions
source("../scripts/pantry.R")


# Load 1 dataInput file & remove food-specific vars 
dat <- readRDS(paste0("../data/processed/rediMR/", pheno_tag, "_datInput.rda"))


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

## diet traits
foods <- c("raw_veg", "cooked_veg", "fresh_fruit", "dried_fruit", "procmeat", "coffee")
names(foods) <- c("Raw vegetables", "Cooked vegetables", "Fresh fruit", "Dried fruit",  "Processed meat", "Coffee")
foods.l <- as.list(foods)
foods.f <- factor(foods, levels=c("raw_veg", "cooked_veg","fresh_fruit", "dried_fruit", "procmeat", "coffee"))


# ===================
## Covariates
# ===================

# covariates to adjust for in ReDiMR
adjCovars <- c("smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab", 
               "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10))

## Format covariates for adjustment 
adjCovarNames <- c(
  smoke_level.lab="Smoking", alch_freq.lab="Alcohol", 
  pa_met_excess_level.lab="Physical Activity", income_level.lab = "Income", 
  educ_level.lab="Education", bmi="BMI", waist2hip="Waist-to-hip",
  dietPC1="Diet Pattern PC1", dietPC2="Diet Pattern PC2", dietPC3="Diet Pattern PC3",
  dietPC4="Diet Pattern PC4", dietPC5="Diet Pattern PC5", 
  #`dietPC1+dietPC2+dietPC3+dietPC4+diePC5`="Top 5 Diet Patterns", 
  dietPC6="Diet Pattern PC6", dietPC7="Diet Pattern PC7", dietPC8="Diet Pattern PC8", 
  dietPC9="Diet Pattern PC9", dietPC10="Diet Pattern PC10"
)



#################################################################################
##  Descriptive plots of diet traits, diet patterns & covariates 
################################################################################


## Panel FOOD Histograms
plots_hist <- dat %>% select(IID, all_of(foods)) %>%
  pivot_longer(-IID, names_to="Food") %>% 
  mutate(Food=factor(Food, levels=names(foods))) %>%
  ggplot(aes(x=value)) + ggtheme + theme_bw() + 
  facet_grid(~Food, scales = "free") + 
  geom_histogram(bins=10, alpha=0.75) + 
  xlab("Food intake, tablespoons/day") + ylab("Frequency") 

# --> conbine & save as ggplot
ggsave(filename=paste0(out_dir, "/foods_plotHist.pdf"), plots_hist, height=5, width=10)



# =========================================
##  Rotated factor loadings for diet PCs
# =========================================

dietPC_rotations.l <- lapply(foods.l, function(food) {
  (readRDS(paste0("../data/processed/rediMR/", food , "_v2_dietPCs.rda")))$rotation
  #rbind(d, matrix(NA,1,ncol(d), dim=list(food,NULL)))
  }) ; names(dietPC_rotations.l) <- foods

## waterfall plot of Factor loading
plot_dietPC.fun <- function(food, nPCs=10) {
  as.data.frame(dietPC_rotations.l[[food]]) %>%
    mutate(Diet=diet_labels[gsub("_QT", "", gsub("_BIN", "", rownames(.)))]) %>%
    select(Diet, c(paste0("PC", 1:nPCs))) %>%
    mutate(Diet=factor(Diet, levels=Diet[order(PC1, decreasing=T)])) %>%
    pivot_longer(-Diet) %>% 
    mutate(name=factor(name, levels=c(paste0("PC", 1:nPCs))),
           dir_col = factor(case_when(
             value <= -0.2~"Negative/Major", value > -0.2 & value <0 ~ "Negative/Minor",
             value > 0 & value <0.2 ~"Positive/Minor", value >= 0.2 ~ "Positive/Major"),
             levels=c("Negative/Major", "Negative/Minor", "Positive/Minor", "Positive/Major"))) %>%
    ggplot(aes(x=value, y=Diet, fill=dir_col)) + ggtheme + 
    facet_grid(~name) + geom_col() + ylab(" ") + 
    xlab(paste("Rotated Factor Loadings")) +
    geom_vline(xintercept = c(-0.2, 0.2), color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values=palettes$pugr, name = "Factor Loading") +
    theme(legend.position = "bottom") +
    ggtitle(paste("PCA-deried dietary patterns without: ", diet_labels[food])) 
}

for(food in foods) {
  ggsave(filename = paste0(out_dir,"/", food, "_plotDietPCs.pdf"), plot_dietPC(food), height=4, width=12)
}

for(food in foods){
  ggsave(plot_dietPC.fun(food,5), filename=paste0("../output/",food,"_dietPCs_top5.pdf"), height=4, width=8)
}


# =========================================
## FOOD x covariate plots 
# =========================================

cat("Generating descriptive plots of covariates for adjustment. \n ")

# compile covariates
catCovars <- c(adjCovarNames[c("smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab",
                           "income_level.lab", "educ_level.lab")], 
               "bmi_tert"="BMI", "waist2hip_tert"="Waist-to-hip")

# Functions 
plot_dietXcatcov <- function(food, catCovariates) {
  
  dat <- readRDS(paste0("../data/processed/rediMR/", food, "_v2_datInput.rda")) %>% mutate(
    bmi_tert = cut(bmi, breaks=quantile(bmi, probs=seq(0,1,0.33), na.rm=T), 
                   labels=c("Tertile_1", "Tertile_2", "Tertile_3")),
    waist2hip_tert = cut(waist2hip, breaks=quantile(waist2hip, probs=seq(0,1,0.33), na.rm=T), 
                         labels=c("Tertile_1", "Tertile_2", "Tertile_3")))
  
  barDat <- dat %>% select(y=all_of(food[[1]]), all_of(names(catCovars))) %>% 
    filter(complete.cases(y)) 
  do.call(rbind.data.frame, lapply(names(catCovars), function(cov) {
    d <- as.data.table(do.call(rbind, tapply(barDat$y, barDat[cov], function(x) cbind(m=mean(x), sd=sd(x)))))
    d <- d %>% mutate(Covariate=rep(catCovars[cov], nrow(d)), Levels=names(table(barDat[cov])),
                    Levels=factor(Levels, levels=c(Levels[1:nrow(.)])), Col=factor(1:nrow(.)))
    return(d) })) %>%
    mutate(Covariate=factor(Covariate, levels=c("Smoking", "Alcohol", "Physical Activity", "Income", "Education", "BMI", "Waist-to-hip"))) %>%
    ggplot(aes(x=Levels, y=m, ymin=m, ymax=m+sd, fill=Covariate)) + 
    facet_grid(~Covariate, scale="free", space="free_x") + ggtheme + 
    theme(axis.text.x=element_text(angle=35, hjust=1, color="black")) +
    geom_bar(stat="identity", width=0.75) + geom_errorbar(width=0.15) +
    scale_fill_manual(values=palettes$CovarGroups) +
    theme(legend.position = "none", axis.text.x = element_text(angle=35, size=7, color="black") ) +
    xlab(" ") + ylab("Mean (SD) Intake") +
    ggtitle(paste("Mean (SD)", diet_labels[food], "by Lifestyle & SES Covariates"))
}


################################################################################
##  Build rediMR plots of %B change for ALL & EACH covariate adjustments
################################################################################

# =====================================
## Adjusting for ALL covariates  
# =====================================

cat (paste0(" \n Making ReDiMR plots. \n "))

# Dot plot of pctBchange when adjusting for ALL covariates
plot_BChange_all.fun<-function(food, tag) {
  tab_all <- fread(paste0(rediMR_dir, "/", food, "_", tag, "_tabBchangeAllCov.csv")) %>% 
    arrange(B_pctChange) %>%
    mutate(snp=gsub("snp","",gsub("[.]",":",snp))) %>%
    mutate(SNP=factor(snp, levels=snp), RefinedSet=factor(RefinedSet, levels=c("0","1")))
  xscale<-max(abs(ceiling(tab_all$B_pctChange*1.10)))
    
  tab_all %>% 
    ggplot(aes(x=B_pctChange, y = SNP, color = RefinedSet)) + 
    theme_bw() + ggtheme + 
    scale_x_continuous(limits=c(-xscale, xscale)) +
    geom_vline(xintercept = 0, color = "black") + 
    geom_vline(xintercept = c(-20, 20), color = "black", linetype = "dashed") + 
    geom_point(size=2.5) + 
    scale_color_manual(values = c(palettes$NatComms[1], palettes$NatComms[3]), 
                       name = "Refined Set", labels=c("Excluded", "Included")) +
    xlab("Beta magnitude (%) change from Base model") + ylab(" ") +
      ggtitle(paste0("ReDiMR Plot: B pct change for all covariates")) + 
      theme(title = element_text(size=8), legend.position = "bottom",
            plot.title = element_text(hjust=-0.05),
            axis.text.y=element_text(size=6)) 
}

plot_BChange_all.fun("fresh_fruit","v2")

# ================================
## Adjusting for EACH covariate
# ================================

# Dot plot of pctBchange when adjusting for EACH covariate

plot_Bchange_each.fun <- function(food, tag) {
  tab_each <- fread(paste0(rediMR_dir, "/", food, "_", tag, "_tabBchangeByCov.csv")) %>% 
    mutate(covar=factor(covar, levels=adjCovarNames)) %>%
    arrange(covar) 
  yscale <- ceiling(max(abs(tab_each$B_pctChange))*1.10)
  tab_each %>%
    ggplot(aes(x=covar, y = B_pctChange, color = covar)) + 
    theme_bw() + ggthemeF + 
    theme(axis.text.x = element_text(angle=35, vjust=0.99, size=8),
          legend.position = "none") +
    #scale_y_continuous(limits=c(-yscale, yscale)) +
    geom_hline(yintercept = 0, color = "black") + 
    geom_hline(yintercept = c(-10, 10), color = "black", linetype = "dashed") + 
    geom_point(size=1.15, position = position_jitter(0.25)) + 
    scale_color_manual(values = palettes$CovarGroups, name = "Adjusted\nCovariate") +
    #scale_y_discrete(labels=rev(plot_B_pct_change_full_model_dat$snp)) +
    ylab("Beta magnitude (%) change from Base model") + xlab(" ") + 
    ggtitle(paste0("ReDiMR Plot: B pct change for each covariate")) +
      theme(title = element_text(size=8),
            axis.title.y = element_text(color="black", size=8))
}



################################################################################
## Compile panel summary figure
################################################################################

# Function to compile panel figure by food
plot_redimr_summary.fun <- function(food) {
  ggarrange(ggarrange(plot_dietXcatcov(food, catCovars),
                      plot_Bchange_each.fun(food, "v2"), ncol=1, heights=c(0.75, 1.25),labels=c("A. ", "B. ")),
            plot_BChange_all.fun(food, "v2"), labels=c("", "C."), ncol=2, widths=c(1.5, 1))
}

plot_redimr_summary.fun(foods[6])

# plot descriptives
for(food in foods) {
  ggsave(plot_dietXcatcov(food), filename=paste0("../output/",food,"_v2_plotCovars.pdf"), height=3, width=12)
}


# plot Bchange by covariates
for(food in foods) {
  ggsave(plot_Bchange_each.fun(food, "v2"), filename=paste0("../output/",food,"_v2_plotBchangeByCov.pdf"), height=4, width=8)
}


# plot Bchange all covariates
for(food in foods) {
  ggsave(plot_BChange_all.fun(food, "v2"), filename=paste0("../output/",food,"_v2_plotBchangeAllCov.pdf"), height=8, width=5)
}

#EOF

ggsave(plot_BChange_all.fun("raw_veg", "v2"), filename=paste0("../output/","raw_veg","_v2_plotBchangeAllCov.pdf"), height=5, width=5)
ggsave(plot_BChange_all.fun("procmeat", "v2"), filename=paste0("../output/","procmeat","_v2_plotBchangeAllCov.pdf"), height=2.5, width=5)

