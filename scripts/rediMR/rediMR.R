# ReDiMR - v1
# Last updated: May 15, 2023



############
## Set Up ##
############

# load required packages
library(tidyverse) ; library(data.table) ; library(dplyr) ; library(parallel)
library(paletteer) ; library(RColorBrewer)

# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
tag <- args[2]

pheno.tag <- paste0(pheno, ".", tag)
ssInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno.tag, "_ssInput.csv")  #args[2]
datInput <- paste0("../data/processed/rediMR/", pheno, "/", pheno.tag, "_datInput.rda") #args[3]
outDir <- dirname(datInput) #args[4]
pctBdeltIncl <- 20 #args[5]


# load basic functions
source("../scripts/basic_functions.R")

# create directory to store results
system(paste0("mkdir -p ", outDir))



## default ReDiMR variable assignments // optional arguments ?

# covariates in base gwas
gwasCovars <- c("age","sex", paste0("gPC", 1:10))

# covariates to adjust for in ReDiMR
adjCovars <- c("smoke_level.lab", "alch_freq.lab", "pa_met_excess_level.lab", 
               "income_level.lab", "educ_level.lab", "bmi", "waist2hip", paste0("dietPC", 1:10))



######################################
## Format covariates for adjustment ##
######################################

# Assign covar names
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

# Write as lists
adjCovars.l <- as.list(adjCovars)
adjCovarNames.l <- as.list(adjCovarNames)

nCovars.l <- as.list(1:length(adjCovarNames.l))



############################################################
## Build function to calculate pctBchange from Base Model ##
############################################################

pctBchange.fun <- function(pheno, snp, adjCovar, baseCovars=gwasCovars, replace_covar_name=NA, data=dat) {
  
  covarsFormatted <- paste0(baseCovars, collapse="+")
  snp_EA <- names(data %>% select(starts_with(snp)))
  baseM <- lm(formula(paste0(pheno, "~", snp_EA, "+", covarsFormatted)), data)
  adjM <- lm(formula(paste0(pheno, "~", snp_EA, "+", covarsFormatted, "+", adjCovar)), data)
  
  baseCI <- confint(baseM)[2,]
  adjCI <- confint(adjM)[2,]
  pctB <-round(((adjM$coef[2]-baseM$coef[2])/baseM$coef[2])*100, 2)
  
  baseP=summary(baseM)$coef[2,4]
  adjP=summary(adjM)$coef[2,4]

  if(!is.na(replace_covar_name)) {covarName = replace_covar_name} else {covarName = adjCovar}
  lmOut <- cbind.data.frame(
    snp=paste(snp), covar=covarName, B_base=baseM$coef[2], B_adj=adjM$coef[2], 
    lowCI_base=baseCI[1], upCI_base=baseCI[2], lowCI_adj=adjCI[1], upCI_adj=adjCI[2], 
    P_base=baseP, P_adj=adjP, B_pctChange=pctB) ; rownames(lmOut) <- paste0(covarName)
  
  return(lmOut)
}



###############################################
## Load input files & Run basic descriptives ##
###############################################

# Load summary stats data
ss<-fread(ssInput) %>% filter(LOCI == 1)

# Load phenotype/dosage data
dat<-readRDS(datInput)

# list of snps
snps <- names(dat %>% select(contains(ss$SNP)))



#######################################
## Basic descriptives by pheno_quant ##
#######################################

tab_descrByPhenoQ <- do.call(rbind.data.frame, lapply(nCovars.l, function(i) {
  tmp <- dat %>% select(phenoQ, cov=all_of(names(adjCovarNames)[i])) %>%
    filter(!is.na(phenoQ)) %>% group_by(phenoQ)
  if(is.factor(tmp$cov)) {
    out <- do.call(rbind.data.frame, lapply(levels(tmp$cov), function(lvl) {
      tmp %>% summarise(cov_level = n_pct(cov, level = lvl)) %>% 
        t() %>% as.data.frame() %>% filter(V1 != "Q1") } ))
    rownames(out) <- paste0( adjCovarNames[i], "_", levels(tmp$cov), ", n (%)") } 
  else if(is.numeric(tmp$cov)) {
    out <- tmp %>% summarise(cov_msd = mean_sd(cov)) %>%
      t() %>% as.data.frame() %>% filter(V1 != "Q1")
    rownames(out) <- paste0(adjCovarNames[i], ", mean \u00B1 SD") }
  colnames(out) <- c(paste0("Q", 1:ncol(out)))
  out
} )) ; head(tab_descrByPhenoQ)


cat(paste0("Writing file with descriptive characteristics by ", pheno, " quantiles: \n"))
print(dat %>% group_by(phenoQ) %>% select(phenoQ, pheno=all_of(pheno)) %>% summarise(m_SD=mean_sd(pheno, d=3)))

# Write descriptives to csv
fwrite(tab_descrByPhenoQ, file = paste0(outDir, "/", pheno, "_descrByPhenoQs.csv"))



############################################################
## Calculate pct % beta for ALL covariates --> Refinement ##
############################################################

cat("Calculating pctBchange when adjusting for ALL covariates ... ")

# For eacn SNP, tabulate pctBchange when adjusting for ALL covariates
tab_pctBchangeAllCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) {
  pctBchange.fun(pheno=pheno, snp=snp, adjCovar=paste0(adjCovars, collapse="+"), 
                 replace_covar_name = "All_Covariates", data=dat)}, mc.cores = 8 )) %>% 
  mutate(RefinedSet = ifelse(abs(B_pctChange) < pctBdeltIncl,1,0)) 


# Write results to csv
fwrite(tab_pctBchangeAllCov, file = paste0(outDir, "/", pheno.tag, "_pctBchangeAllCov.csv"))

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno.tag, "_pctBchangeAllCov.csv")))
head(tab_pctBchangeAllCov)



###############################################################
## Calculate pct % beta BY covariate --> Covariate Influence ##
###############################################################

cat("Calculating pctBchange when adjusting for EACH covariate ... ")

# For each snp, tabulate pctBchange when adjusting for EACH covariate
tab_pctBchangeByCov <- do.call(rbind.data.frame, mclapply(snps, function(snp) { 
  do.call(rbind.data.frame, mclapply(nCovars.l, function(i) {
    pctBchange.fun(pheno=pheno, snp=snp, adjCovar=names(adjCovarNames)[i], 
                   replace_covar_name=adjCovarNames[[i]], data=dat) }, mc.cores = 8))  }, mc.cores = 8))


# Write results to csv
fwrite(tab_pctBchangeByCov, file = paste0(outDir, "/", pheno.tag, "_pctBchangeByCov.csv"))

cat (paste0("DONE: Results written to ", paste0(outDir, "/", pheno.tag, "_pctBchangeByCov.csv")))
head(tab_pctBchangeByCov)




###################
## ReDi MR Plots ##
###################

cat (paste0("Making ReDiMR plots."))

palettes <- list(NatComms=paletteer_d("ggsci::nrc_npg", n=10))

# set default ggplot theme
ggthemeF <- theme(panel.grid.minor.y = element_blank(), 
                  panel.grid.minor.x = element_blank(), 
                  axis.text = element_text(color="black", hjust=1),
                  #axis.title = element_text(face = "bold", size=11, vjust=1.5),
                  #plot.title=element_text(size=12),
                  #legend.text = element_text(size=10), legend.title = element_text(face="bold", size=10),
                  legend.position = "right", legend.box.background = element_rect(color = "black"))


# Dot plot of pctBchange when adjusting for ALL covariates
xscale <- ceiling(max(abs(tab_pctBchangeAllCov$B_pctChange))*1.15)

dotAll <- tab_pctBchangeAllCov %>%
  arrange(B_pctChange) %>%
  mutate(SNP=factor(snp, levels=snp), RefinedSet=factor(RefinedSet, levels=c("0","1"))) %>%
  ggplot(aes(x=B_pctChange, y = SNP, color = RefinedSet)) + 
  theme_bw() + ggthemeF + 
  scale_x_continuous(limits=c(-xscale, xscale)) +
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = c(-20, 20), color = "black", linetype = "dashed") + 
  geom_point(size=2.5, position = position_jitter(0)) + 
  scale_color_manual(values = c(palettes$NatComms[1], palettes$NatComms[3]), 
                     name = "Refined Set", labels=c("Excluded", "Included")) +
  #scale_y_discrete(labels=rev(plot_B_pct_change_full_model_dat$snp)) +
  xlab(" % |Beta| change from Base model") + ylab(" ") +
  ggtitle(paste0("Refined set of loci for ", pheno, "\n<", pctBdeltIncl, "% change in B"))


# Dot plot of pctBchange when adjusting for EACH covariate

palettes$CovarGroups=c(brewer.pal(9,"Oranges")[4:6], brewer.pal(9, "Greens")[5:4], 
   brewer.pal(9, "Blues")[5:6],  brewer.pal(9, "Purples")[c(8:4)], brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6])

yscale <- ceiling(max(abs(tab_pctBchangeByCov$B_pctChange))*1.10)
dotCov <- tab_pctBchangeByCov %>%
  mutate(covar=factor(covar, levels=adjCovarNames)) %>%
  arrange(covar) %>%
  ggplot(aes(x=covar, y = B_pctChange, color = covar)) + 
  theme_bw() + ggthemeF + 
  theme(axis.text.x = element_text(angle=35, vjust=0.99, size=10),
        legend.position = "none") +
  scale_y_continuous(limits=c(-yscale, yscale)) +
  geom_hline(yintercept = 0, color = "black") + 
  geom_hline(yintercept = c(-10, 10), color = "black", linetype = "dashed") + 
  geom_point(size=1.15, position = position_jitter(0.25)) + 
  scale_color_manual(values = palettes$CovarGroups, name = "Adjusted\nCovariate") +
  #scale_y_discrete(labels=rev(plot_B_pct_change_full_model_dat$snp)) +
  ylab(" % |Beta| change from Base model") + xlab(" ") + 
  ggtitle(paste0(pheno, " % |Beta| change by SNP and covariate"))



# Save plots as PDF
pdf(paste0(outDir, "/", pheno.tag, "_plotBpctAllCov.pdf"), height = 5.5, width = 6)
dotAll
dev.off()

pdf(paste0(outDir, "/", pheno.tag, "_plotBpctByCov.pdf"), height = 4.5, width = 6)
dotCov
dev.off()

cat(paste0("Done making plots written to ", outDir, "/", pheno.tag, "_plotOutputs.pdf"))

cat("
Congratulatulations!! You completed ReDiMR Step 1 (SNP Refinement).

    --> You are now ready to proceed to ReDiMR Step 2-Mendelian Randomization. Good luck!")



## EOF - TEMPORARY

#*add quadrant plot??
#*

