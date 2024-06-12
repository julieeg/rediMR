#"RediMR plots for testsnps 




################################
##  Set up & Load data files  ##
################################

# =======================
##  Set up
# =======================

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno=args[1]
tag="testsnps"
pheno_tag=paste0(pheno,"_",tag)

testDir="../data/processed/testsnps"
datInput=paste0(testDir, "/", pheno_tag, "_datInput.rda")

# load basic functions
source("../scripts/basic_functions.R")


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


# ========================
## Test SNPs
# ========================

test_snps <- c("rs6690619_C", "rs7619139_T", "rs1421085_T", "rs838133_A",  "rs9323534_C",  "rs1726866_G")
names(test_snps) <- c("GATAD2B_SES", "RARB_GeneralDiet", "FTO_Adiposity", "FGF21_BiolORSweet", "OR4K17_SmellPercept", "TAS2R38_TastePercept")

test_snps_labs <- c(rs6690619_C="rs6690619_C - GATAD2B_SES", rs7619139_T="rs7619139_T - RARB_GeneralDiet", 
                    rs1421085_T="rs1421085_T - FTO_Adiposity",  rs838133_A="rs838133_A - FGF21_BiolORSweet", 
                    rs9323534_C="rs9323534_C - OR4K17_SmellPercept",  rs1726866_G="rs1726866_G - TAS2R38_TastePercept")



#################################################################################
##  Descriptive plots of diet traits, diet patterns & covariates 
################################################################################

# set ggplot theme
ggtheme <- theme(panel.grid.minor.y = element_blank(), 
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




# =========================================
##  Rotated factor loadings for diet PCs
# =========================================

dietPCs <- readRDS(paste0("../data/processed/testsnps/",pheno_tag,"_dietPCs.rda"))$rotation

## waterfall plot of Factor loading
plot_dietPC.fun <- function(nPCs=10) {
  as.data.frame(dietPCs) %>%
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
    theme(legend.position = "bottom")
}

plot_dietPC.fun(5)
ggsave(plot_dietPC.fun(10), filename = paste0(testDir,"/", pheno_tag, "_plotDietPCs.pdf"), height=4, width=12)





# ===========
## % beta change by covariate
# =========

tabBchangeAllCov <- fread(paste0(testDir, "/",pheno_tag,"_tabBchangeAllCov.csv"))
tabBchangeByCov <- fread(paste0(testDir, "/",pheno_tag,"_tabBchangeByCov.csv"))

tab_each <- tabBchangeByCov %>%
  mutate(covar=factor(covar, levels=adjCovarNames)) %>%
  arrange(covar)
#yscale <- ceiling(max(abs(tab_each$B_pctChange))*1.10)


tab_each %>%
  #filter(covar != "Diet Pattern PC1") %>%
  ggplot(aes(y=covar, x=B_pctChange, color=covar, group=covar)) + 
  theme_bw() + ggthemeF + 
  facet_wrap(~snp, scales="free", labeller = as_labeller(test_snps_labs)) +
  xlab("Covariates") + ylab("% change in |Beta| from Base model") +
  geom_point(position = position_dodge(0.15), size=3.5) + #position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1)
  geom_vline(xintercept = 0, color = "black", linewidth=0.5) +
  scale_color_manual(values = palettes$CovarGroups, name = "Adjusted\nCovariate") +
  geom_vline(xintercept = seq(-30,30,20), color="#00000015") +
  geom_vline(xintercept = c(-10,10), color="#000000", linetype = "dashed")


