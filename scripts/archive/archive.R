# Archive of retired plotting functions




# ==========================================================================
## Function to plot BETA (adjusted) when adjusting for each CONFOUNDER
# ==========================================================================

plot_adjB_conf.fun <- function(Bchange_conf.l) {
  
  tab_conf <- Bchange_conf.l$BchangeBy %>%
    pivot_longer(c(ends_with("base"), ends_with("adj")), names_to=c("stat", "mod"), names_sep="_") %>%
    pivot_wider(names_from="stat", values_from = value) %>%
    mutate(Covar = ifelse(Covar == "Smoking" & mod == "base", "Unadjusted", Covar)) %>%
    filter(Covar == "Unadjusted" | mod == "adj") %>% 
    mutate(Covar = factor(Covar, levels=c("Unadjusted", covarSets$stnd$Names))) %>%
    mutate(B_pctChange = ifelse(Covar == "Unadjusted", 0, B_pctChange)) 
  
  yscale<-max(abs(tab_conf$B))*1.05
  
  tab_conf %>%
    ggplot(aes(x=snp, y = B, color = Covar)) + 
    theme_bw() + ggtheme2 +
    theme(axis.text.x = element_text(angle=35, hjust=1, size=8),
          panel.grid.major.x = element_line()) +
    geom_point(size=1.5, position=position_jitter(0.35), alpha=0.85) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
    scale_y_continuous(limits=c(-yscale, yscale)) +
    scale_color_manual(values =  c("black",palettes$Conf), name = "Confounder") +
    ylab("Beta (95% CI)") + xlab(" ")
}

plot_adjB_conf.fun(veg_Bchange_conf.l)
plot_adjB_conf.fun(fruit_Bchange_conf.l)



##### Assumption 1b: Confounder-diet associations become attenuated after adjusting for diet PCs

## Function to plot trait-confounder associations
plot_CorrConf.fun <- function(pheno, tab) {
  xlim=c(0, length(unique(tab$covar)))
  tab %>% 
    pivot_longer(c("B")) %>%
    mutate(nPCs = as.numeric(gsub("Top", "", gsub("DietPCs.*", "", gsub("Unadjusted", "0", covar))))) %>%
    mutate(statName=factor(name, levels=c("B", "adjR2", "Fstat"), 
                           labels = c("Beta (95% CI)", "Model Adjusted R2", 
                                      "Model F-statistic"))) %>%
    mutate(lowCI=ifelse(name=="B", lowCI, NA),
           upCI=ifelse(name=="B", upCI, NA),
           col=ifelse(nPCs==0, "A", "B")) %>%
    mutate(conf=factor(conf, levels=confounder_vars.num)) %>%
    ggplot(aes(x=nPCs, y=value, ymin=lowCI, ymax=upCI, shape=statName, 
               color=conf, fill=conf)) + 
    facet_wrap(~conf, scales="free_y", nrow=2, labeller = as_labeller(confounders_num.l)) +
    #facet_wrap(~statName, scales="free_y", ncol=2) +
    geom_hline(yintercept = 0, linetype="dashed", linewidth=0.15) +
    #geom_vline(xintercept = c(1,5,10,15,20), linetype="solid", linewidth=0.1) +
    geom_point(size=1.5, alpha=85) + #, color="grey38"
    geom_errorbar(width=0.25, linewidth=0.25) +
    scale_x_continuous(breaks=seq(0,23,3), limits=c(xlim[1], xlim[2])) + 
    geom_line() +
    scale_shape_manual(values=c(15, 17, 19), guide=F) +
    scale_color_manual(values=palettes$StndConf_num, name="") +
    scale_fill_manual(values=palettes$StndConf_num, name="") +
    ylab("Beta (95% CI)") + xlab("") +
    ggtheme +
    theme(panel.grid = element_blank(), panel.grid.minor.y = element_blank(),
          legend.position = "none", 
          axis.title.y = element_text(face="bold"),
          axis.text.x = element_text(size=8))
}


# ==========================
## Raw vegetable intake 
# ==========================

# dietGroup2 (dPCs WITHOUT the diet trait)
plot_CorrConf.fun("raw_veg", veg_confPheno$dietGroup2) #%>%
#ggsave(filename="../output/p_veg_v2_dietGroup2_confPheno_.pdf", height=4, width=9)

# dietGroup3 (dPCs WITH the diet trait)
#plot_CorrConf.fun("raw_veg", veg_confPheno$dietGroup3) %>%
#  ggsave(filename="../output/p_veg_v2_dietGroup3_confPheno.pdf", height=4, width=9)


# ==========================
## Fresh fruit intake 
# ==========================

# dietGroup2 (dPCs WITHOUT the diet trait)
plot_CorrConf.fun("fresh_fruit", fruit_confPheno$dietGroup2) #%>%
#ggsave(filename="../output/p_fruit_v2_dietGroup2_confPheno_.pdf", height=4, width=9)

# dietGroup3 (dPCs WITH the diet trait)
#plot_CorrConf.fun("fresh_fruit", fruit_confPheno$dietGroup3) %>%
#  ggsave(filename="../output/p_fruit_v2_dietGroup3_confPheno.pdf", height=4, width=9)






## % change in beta for EACH dietPC =============
plot_Bchange_dPC <- plot_Bchange_dPCconf.fun(veg_Bdat.l$each_dPC, "raw_veg", varset="dPCs") +
  scale_x_discrete(labels=seq(1,23,1)) + 
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) + 
  xlab("# dietPCs adjusted")
#plot_Bchange_dPC %>% ggsave(filename = "../output/plot_raw_veg_v2_Bchange_dPC.pdf", height = 4, width=6)


## adjB for each confounder  =============
tab_veg_covarEach_adjB <- veg_Bdat.l$each_conf %>%
  pivot_Bdat_to_long() %>% 
  mutate(shape=ifelse(Covar=="Unadjusted", "square", "dot")) %>%
  mutate(Covar=factor(Covar, levels=c("Unadjusted", covarSets$stnd$Names)))

plot_adjB_conf <- tab_veg_covarEach_adjB %>%
  ggplot(aes(x=Covar, y = B, group=snp, color = Covar, shape=shape)) + 
  theme_bw() + ggtheme2 + 
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_point(size=1.55, alpha=0.95) +
  geom_line(linewidth=0.15) +
  scale_color_manual(values=c("#00000065", palettes$Conf)) +
  scale_shape_manual(values=c("square"=15, "dot"=19)) +
  ylab("Beta (95% CI)") + 
  #xlab("Covariates Adjusted") +
  xlab(NULL) +
  theme(legend.position = "none") 

plot_adjB_conf




# ==============================================================================
## Dot plot BETA (% change) when adjusting for EACH covariate ON TOP of dietPCs
### Note: dietPCs PLUS confounders
# ==============================================================================

tab <- rbind.data.frame(
  veg_Bdat.l$each_conf, veg_Bdat.l$each_dPCplusconf, 
  veg_Bdat.l$all_conf %>% mutate(Covar=gsub("Standard", "All Confounders", Covar)),
  veg_Bdat.l$all_dPCplusconf %>% filter(Covar =="DietPCs_Add7cWaist2Hip")) %>% 
  left_join(veg_Bdat.l$all_dPCs %>% select(c("snp", B_dPC="B_adj")), by = "snp") %>%
  rowwise() %>%
  mutate(B_pctChange = ifelse(startsWith(Covar, "DietPCs_"), ((B_adj-B_dPC)/B_dPC)*100, B_pctChange),
         Reference_Model = ifelse(startsWith(Covar, "DietPCs_"), "DietPC", "Base")) %>%
  pivot_Bdat_to_long() %>% 
  mutate(Covar = factor(Covar, levels=c("Unadjusted", "DietPCs", "Smoking", "DietPCs_Smoke", "Alcohol", "DietPCs_Alch", "Physical Activity", "DietPCs_PhysAct",
                                        "Income", "DietPCs_Income", "Education", "DietPCs_Educ", "BMI", "DietPCs_BMI", "Waist-to-hip", "DietPCs_Waist2Hip",
                                        "All Confounders", "DietPCs_Add7cWaist2Hip"),
                        labels=c("Unadjusted", "DietPCs", "Smoking", "DietPCs+Smoking", "Alcohol", "DietPCs+Alch", "Phys Act", "DietPCs+Phys Act",
                                 "Income", "DietPCs+Income", "Education", "DietPCs+Education", "BMI", "DietPCs+BMI", "Waist2Hip", "DietPCs+Waist2Hip",
                                 "All Confounders", "DietPCs+All Confounders"))) %>%
  mutate(Covar_cat = gsub("DietPCs[+]", "", Covar)) %>%
  mutate(Covar_cat = ifelse(Covar_cat == "Alch", "Alcohol",ifelse(Covar_cat == "Add7cWaist2Hip", "All Confounders", Covar_cat))) %>%
  mutate(Covar_cat = factor(Covar_cat, levels=c("Unadjusted", c("Smoking", "Alcohol", "Phys Act", "Income", "Education", "BMI", "Waist2Hip", "All Confounders"))))


yscale <- ceiling(max(abs(tab$B_pctChange))*1.02)

plot_Bchange_confCompare <- tab %>%
  filter(Covar != "Unadjusted") %>%
  ggplot(aes(x=Reference_Model, y = B_pctChange, color = Covar, group=snp, shape=Reference_Model)) + 
  theme_bw() + ggtheme2 +
  facet_grid(~Covar_cat, scales = "free_x") +
  # theme(axis.text.x = element_text(angle=35, hjust=1)) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_point(size=1.55, alpha=0.95) + 
  geom_line(linewidth=0.15) +
  geom_hline(yintercept = c(-10, 10), color = "black", linetype = "dashed") + 
  scale_y_continuous(limits = c(-yscale, yscale)) +
  scale_color_manual(values =  c(rep(palettes$Conf[1:7], each=2), "grey38", "grey38"), guide=F) +
  ylab("% change in Beta") + xlab(" ") + 
  theme(legend.position = "none")

plot_Bchange_confCompare %>%
  ggsave(filename="../output/raw_veg_v2_plot_panel_Fig2_Bchange_compare.pdf", height = 3, width =8)

plot_Bchange_confCompare

