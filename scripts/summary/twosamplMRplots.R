#twosampleMRplots.R


#################################
## Set up & assign parameters  ##
#################################

# =======================
##  Set up
# =======================

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)

#remotes::install_github("MRCIEU/TwoSampleMR")
#remotes::install_github('MRCIEU/ieugwasr')
library(TwoSampleMR) ; library(ieugwasr)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
tag <- args[2]

pheno_tag <- paste0(pheno, "_", tag)
rediMR_inputPref=paste0("../data/processed/rediMR/", pheno_tag)

mr_dir="../data/processed/mr"
datInput <- paste0(rediMR_inputPref, "_datInput.rda")  #args[4] 
outDir <- "../output"


# load basic functions
source("../scripts/basic_functions.R")


#################################################################
##  Load input files ; Write functions to calculate %B change  ##
#################################################################

# ==========================================
## Load data input files
# ==========================================

# Load phenotype/dosage data
dat<-readRDS(datInput)

# Load MRsummary_rda files
files <- list.files(path=mr_dir, pattern=paste0(pheno_tag,"_MRsummary"))
mr_summary.l <- lapply(files, function(file) readRDS(paste0(mr_dir,"/",file)))
names(mr_summary.l) <- c(gsub(".*summary_","", gsub(".rda","",files)))

# Make vectors of outcomes & snpsets
outcomes <- names(mr_summary.l)
snpsets <- names(mr_summary.l[[1]])



# ==========================================
## Build ggplot theme
# ==========================================

ggtheme <- theme_bw() + theme(panel.grid.minor.y = element_blank(), 
                 panel.grid.minor.x = element_blank(), 
                 axis.text.x = element_text(color="black", size=8), #angle=35, hjust=1, 
                 axis.text.y = element_text(color="black", size=8),
                 axis.title = element_text(color="black", size=10),
                 legend.position = "right", 
                 legend.box.background = element_rect(color = "black"),
                 legend.key.size = unit(0.75, 'line'),
                 legend.margin = margin(0.2,0.2,0.2,0.2, unit="cm"),
                 legend.text = element_text(size=8.5), 
                 legend.title = element_text(face="bold", size=8),
                 plot.title=element_text(size=8),
                 strip.text = element_text(face="bold", size=10)
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


#################################################
##  Basic MR summary plots (from TwoSampleMR)  ##
#################################################

# ==========================================
## MR scatter plot
# ==========================================

plot_mr_scatter.fun <- function(outcome, snpset) {
  dat <- mr_summary.l[[outcome]][[snpset]]$mr_dat
  res <- mr_summary.l[[outcome]][[snpset]]$mr_res
  mr_scatter_plot(res, dat)[[1]] +
    xlab(paste0("SNP effect on ",diet_labels[[pheno]], " (",snpset," SNPs)")) +
    ggtitle(paste("MR Scatter Plot:\n", diet_labels[[pheno]], "&", gsub(" [||].*","", dat$outcome[1]), "-", snpset, "SNPs")) +
    ggtheme
}

plot_mr_scatter_bySNPset.fun <- function(outcome){
  ggarrange(plot_mr_scatter.fun(outcome, "All"), plot_mr_scatter.fun(outcome, "Refined"), 
            plot_mr_scatter.fun(outcome, "Unrefined"), ncol=3, common.legend = T, legend="top")
}

## Make all scatter plots
for(outcome in outcomes){plot_mr_scatter_bySNPset.fun(outcome)}



# ==========================================
## MR leave-one-out forest plot
# ==========================================

# leave one out plot
plot_mr_leaveoneout.fun <- function(outcome, snpset) {
  dat_sinlge <- mr_singlesnp(mr_summary.l[[outcome]][[snpset]]$mr_dat) %>% 
    mutate(SNP=gsub("All - Inverse variance weighted", "All - IVW",SNP))
  mr_forest_plot(dat_sinlge)[[1]] + ggtheme + theme(legend.position = "none") +
    ggtitle(paste("MR Leave-one-out Forest Plot:\n", diet_labels[[pheno]], "&", gsub(" [||].*","", dat_sinlge$outcome[1]), "-", snpset, "SNPs"))
}

plot_mr_leaveoneout_bySNPset.fun <- function(outcome){
  ggarrange(plot_mr_leaveoneout.fun(outcome, "All"), plot_mr_leaveoneout.fun(outcome, "Refined"), 
            plot_mr_leaveoneout.fun(outcome, "Unrefined"), ncol=3)
}

plot_mr_leaveoneout_bySNPset.fun(outcomes[1])



# ==========================================
## MR total forest plot
# ==========================================

## Compile plots for all outcomes
mr_sum_forest <- do.call(rbind.data.frame, lapply(
  as.list(outcomes), function(outcome) { 
    do.call(rbind.data.frame, lapply(snpsets, function(set) {
    mr_summary.l[[outcome]][[set]]$mr_summary %>% 
        #filter(method == "MR Egger" | method == "Inverse variance weighted" | method == "Wald ratio") } )) %>%
        filter(method == "Inverse variance weighted" | method == "Wald ratio") %>% mutate(method=gsub("Inverse variance weighted","IVW",method)) } )) %>%
      mutate(snp_set = factor(snp_set, levels=c("All", "Refined", "Unrefined"))) %>%
      #mutate(method = factor(method, levels=c("MR Egger", "Inverse variance weighted", "Wald ratio"), labels=c("MR Egger", "IVW", "Wald"))) %>%
      #mutate(method = "IVW", "MR Egger") %>%
      mutate(b=ifelse(type=="or", log10(estimate), estimate)) })) %>%
  mutate(outcome=ifelse(id.outcome=="finn-b-I9_STR_EXH_EXNONE", "Ischemic stroke", ifelse(id.outcome=="ebi-a-GCST90018952", "DBP", outcome)))
  
labs<-sapply(c("All","Refined","Unrefined"), function(x) paste0(x, " (N=",mr_sum_forest$nsnp[mr_sum_forest$snp_set==x][1],")"))

plot_mr_forest_bySNPset_allOutcomes <- mr_sum_forest %>% 
  ggplot(aes(x=b, xmin=b-1.96*se, xmax=b+1.96*se, y=snp_set, group=rev(snp_set), color=snp_set, shape=method)) +
  ggtheme + 
  facet_wrap(~outcome, scales="free",ncol=3) +
    geom_vline(xintercept = 0) + ggtheme +
    geom_point(size=2.5, position=position_dodge(0.75)) + 
    geom_errorbarh(position=position_dodge(0.75), height=0.25, lwd=0.5) +
    scale_color_manual(values=c( palettes$NatComms[4], palettes$NatComms[3], palettes$NatComms[1]),
                       name = "Loci Category", labels=labs) +
    scale_shape_manual(values=c(19,17,16), name = "MR Method") +
    xlab("MR Effect Estimate (95% CI)") + ylab("") +
    theme(axis.text.y = element_text(color = "black", size = 10), 
          legend.position = "bottom") #+
  #ggtitle(paste("MR Forest Plot:", diet_labels[pheno], "wtih all MR outcomes"))

ggsave(plot_mr_forest_bySNPset_allOutcomes, filename=paste0(mr_dir,"/",pheno_tag,"_plotMR_forestBySNPset_allOutcomes.pdf"),
       height=7, width=6)



# ==========================================
## Panel plot of heterogeneity nd pleiotropy stats
# ==========================================

## Heterogeneity 
mr_val <- do.call(rbind.data.frame, lapply(
  as.list(outcomes), function(outcome) { 
    do.call(rbind.data.frame, lapply(snpsets, function(set) {
      mr_summary.l[[outcome]][[set]]$mr_summary %>% 
        select(snp_set, outcome, id.outcome, method, nsnp, starts_with("Q"), starts_with("egger")) %>%
        mutate(egger_lci=egger_intercept-1.96*egger_se, 
               egger_uci=egger_intercept+1.96*egger_se) %>%
        filter(method == "MR Egger" | method == "Inverse variance weighted") } )) %>%
      pivot_longer(c(Q, egger_intercept), names_to="measure") %>%
      mutate(Measure=factor(measure, levels=c("Q", "egger_intercept"), labels=c("Heterogeneity", "Horizontal Pleiotropy"))) %>%
      mutate(snp_set = factor(snp_set, levels=c("All", "Refined", "Unrefined"))) %>%
      mutate(method = factor(method, levels=c("MR Egger", "Inverse variance weighted"), labels=c("MR Egger", "IVW"))) 
    })) %>%
  mutate(outcome=ifelse(id.outcome=="finn-b-I9_STR_EXH_EXNONE", "Ischemic stroke", ifelse(id.outcome=="ebi-a-GCST90018952", "DBP", outcome)))
  
labs<-sapply(c("All","Refined","Unrefined"), function(x) paste0(x, " (N=",mr_val$nsnp[mr_val$snp_set==x][1],")"))

plot_mr_het_bySNPset_allOutcomes <- mr_val %>% 
  filter(measure=="Q") %>%
  ggplot(aes(x=method, y=value, fill=snp_set, group=snp_set)) +
  facet_wrap(~outcome, scales="free_y",ncol=4) + ggtheme +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Cochrane's Q") + xlab("") + 
  scale_fill_manual(values=c(palettes$NatComms[6], palettes$NatComms[7], palettes$NatComms[5]),
                    name = "Loci Category", labels = labs) +
  theme(legend.position = "bottom", 
        strip.text = element_text(face="bold", size=8))

plot_mr_pleio_bySNPset_allOutcomes <- mr_val %>%
  filter(measure=="egger_intercept" & method=="MR Egger") %>% 
  ggplot(aes(x=snp_set, y=value, ymin=value-1.96*egger_se, ymax=value+1.96*egger_se, group=snp_set, fill=snp_set)) +
  facet_wrap(~outcome, scales="free_y",ncol=4) + ggtheme +
  geom_bar(stat = "identity", position=position_dodge(0.99)) + 
  geom_hline(yintercept = 0, color = "black", linewidth=0.15) +
  geom_errorbar(width=0.15, position=position_dodge(0.99), linewidth=0.35) +
  ylab("MR Egger Intercept (95% CI)") + xlab("") + 
  scale_fill_manual(values=c(palettes$NatComms[6], palettes$NatComms[7], palettes$NatComms[5]),
                    name = "Loci Category", labels = labs) +
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=25, size=8, hjust=1),
        strip.text = element_text(face="bold", size=8))

ggsave(plot_mr_het_bySNPset_allOutcomes, filename=paste0(mr_dir,"/", pheno_tag,"_plotMR_hetBySNPset_allOutcomes.pdf"), height=3.5, width=6.5)
ggsave(plot_mr_pleio_bySNPset_allOutcomes, filename=paste0(mr_dir,"/", pheno_tag,"_plotMR_pleioBySNPset_allOutcomes.pdf"), height=3.5, width=6.5)



# ==========================================
## Arrange summary figures
# ==========================================

#ggsave(ggarrange(plot_mr_forest_bySNPset_allOutcomes, 
#          ggarrange(plot_mr_het_bySNPset_allOutcomes, plot_mr_pleio_bySNPset_allOutcomes, nrow=2, common.legend = T, legend="none"),
#          ncol=2, widths=c(0.65,1.25), common.legend = T, legend="bottom"),
#       filename=paste0(mr_dir, "/", pheno_tag, "_plotMRpanelBySNPsetAllOutcomes.pdf"), height=6, width=12)



# ==========================================
## Fstat ????
# ==========================================

#cbind(mr_summary.l$colorectal_cancer$All$mr_dat %>% select(SNP, Fstat))
#range(mr_summary.l$colorectal_cancer$Refined$mr_dat %>% select(snp_set, SNP, Fstat))
#mr_summary.l$colorectal_cancer$Unrefined$mr_dat %>% select(snp_set, SNP, Fstat)








