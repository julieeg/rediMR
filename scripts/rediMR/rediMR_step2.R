#Prepare summary stats 


# load required packages
lapply(c("tidyverse", "data.table", "paletteer", "RColorBrewer", "ggpubr", "MendelianRandomization"),  
       library, character.only = TRUE)

# load pantry file with stored parameters, basic functions & ggplot templates
source("../scripts/pantry.R")


#############################################
## ~~~~~~ Set up & input parameters ~~~~~~ ##
#############################################

# ==========================
## Load command arguments
# ==========================
args <- commandArgs(trailingOnly = TRUE)
exposure = args[1] # "oilyfish_QT"
outcome = args[2] # "tg" 
sumstats_file = args[3] # ../data/sumstats/standardized_shared_from_KSJ/std_mrdat_oilyfish_GCST90239664_TG_Graham_GRCh37.csv
covarset = args[4] #"confounders1"
tag = args[5] # "vCole" 
saveDir= args[6] #../data/processed/rediMR/vCole/oilyfish_QT_tg

savePref = paste0(exposure, "_", outcome, "_", covarset)
betachange_file = paste0(saveDir, "/", savePref, "_bchange_full_",tag,".csv")

## Inputs & parameters
inputs <- paste0("\n ReDiMR Step 2 Inputs: \n -exposure ", exposure, 
                 "\n -outcome ", outcome,
                 "\n -sumstats ", sumstats_file,
                 "\n -covarset ", covarset, 
                 "\n -saveDir ", saveDir,
                 "\n -savePref ", savePref,
                 "\n Using files prepared in step 1: ", betachange_file,
                 "\n -tag ", tag, 
                 "\n"
)

cat(paste0(inputs))


# ==========================================
## Load data input files
# ==========================================

## Input files
sumstats <- fread(sumstats_file) %>%
  mutate(Fstat=(beta.exposure^2)/(se.exposure^2))

betachange <- fread(betachange_file) %>%
  mutate_at("snp", ~gsub("id", "", gsub("[.]", ":", .))) %>%
  mutate(snp.exposure = ifelse(startsWith(snp, "rs"), gsub("[_].*", "", snp), sub("(_[^_]*)$", "", snp)))

# Define covariate set
covarSet <- covarSets[[covarset]]

## Merge in rediMR variable for RefinedSets; add flags for differing levels of beta_pctchange
betachange_allcovar <- betachange %>% filter(Covar=="All_Covariates") %>%
  mutate(beta_change.lt20 = ifelse(abs(B_pctChange) < 20, 1, 0),
         beta_change.lt10 = ifelse(abs(B_pctChange) < 10, 1, 0),
         beta_change.lt05 = ifelse(abs(B_pctChange) < 5, 1, 0)
         ) %>%
  select(snp.exposure, Covar, starts_with("beta_change"))

sumstats <- sumstats %>% left_join(betachange_allcovar, by = "snp.exposure")

## Print covariate & outcomes 
outcome_check <- sumstats$outcome[1]
cat("Check outcome in summary stats file vs. input",
    "\n -sumstats: ", outcome_check, 
    "\n -input: ", outcome, "\n"
    )

cat("\nRunning MR for", exposure, "on", outcome, 
    "using harmonized summary statisics (", sumstats_file, "), adjusting for", covarset, ":", 
    paste0("\n   - ", covarSet$Names, " = ", covarSet$Covars), "\n")


## Check number of SNPs in each set and remove any with <2 variants
sumstats %>% reframe(
  "Refined <20%" = n_pct(beta_change.lt20, level="1"),
  "Refined <10%" = n_pct(beta_change.lt10, level="1"),
  "Refined <5%" = n_pct(beta_change.lt05, level="1")
)


# ==========================================
## Run MR for each refinement subset
# ==========================================
mr_dat <- sumstats

## Convert to list of MR datasets
mr_dat.l <- list(
  all=mr_dat,
  refined_lt20=mr_dat %>% filter(beta_change.lt20==1),
  refined_lt10=mr_dat %>% filter(beta_change.lt10==1),
  refined_lt05=mr_dat %>% filter(beta_change.lt05==1)
)

n_snps <- sapply(mr_dat.l, nrow)
if(any(n_snps <= 2)) {
  keep <- which(n_snps > 2)
  mr_dat.l=mr_dat.l[keep]
} else{ keep <- 1:4}
  


#######################################
## Run using Mendelian Randomization ##
#######################################

run_mr.fun <- function(mr_dat) {
  dat <- mr_dat 
  obj <- MendelianRandomization::mr_input(
    bx=dat$beta.exposure, bxse=dat$se.exposure,
    by=dat$beta.outcome, byse=dat$se.outcome,
    exposure="exposure", outcome="outcome", 
    snps=dat$snp.exposure)
  
  res <- MendelianRandomization::mr_allmethods(obj)
  #IVW default is fixed-effects, when 3 or fewer variants
  ivw <- MendelianRandomization::mr_ivw(obj)
  Fstat <- mean(dat$beta.exposure^2/dat$se.exposure^2)
  
  return(list(
    mr_obj=obj,
    mr_res=res,
    mr_ivw=ivw,
    mr_res_df=res@Values, 
    mr_het=data.frame(matrix(
      c(ivw@SNPs,ivw@Heter.Stat, Fstat), 1, 4, dim=list(NULL,c("nSNPs", "Het.Stat", "Het.Stat.Pval", "Fstat"))))
  ))
}

# Run over all subsets of summary stats
mr_output.l <- lapply(mr_dat.l, run_mr.fun)

sets <- names(mr_output.l) 
mr_output.l <- lapply(1:length(mr_output.l), function(set) {
  set <- names(mr_output.l)[set] ; return(
    list(
      mr_obj=mr_output.l[[set]]$mr_obj,
      mr_res=mr_output.l[[set]]$mr_res,
      mr_ivw=mr_output.l[[set]]$mr_ivw,
      mr_res_df=mr_output.l[[set]]$mr_res_df %>% mutate(Set=set),
      mr_het=mr_output.l[[set]]$mr_het %>% mutate(Set=set)
    )) 
}) ; names(mr_output.l)<-sets

for(i in 1:length(mr_output.l)) {
  print(mr_output.l[[i]]$mr_ivw)
}

# =============================================
## Compile & save outputs as data.frames/csvs
# =============================================

set_names=c("All", "Refined <20%", "Refined <10%", "Refined <5%")[keep]
names(sets) <- set_names

## mr_results
mr_output.df <- do.call(rbind.data.frame, lapply(
  mr_output.l, function(set) { set[["mr_res_df"]] } )) %>%
  mutate("exposure" = rep(exposure, nrow(.)),
         "outcome" = rep(outcome, nrow(.)),
         "pair" = rep(paste0(exposure, "_", outcome)),
         "covarset" = rep(covarset, nrow(.)))
colnames(mr_output.df) <- c("Method", "Estimate", "SE", "low_95CI", "up_95CI", "P", "Set", 
                            "exposure", "outcome", "pair", "covarset") 
mr_output.df %>% write.csv(paste0(saveDir, "/", savePref, "_mr_results_", tag, ".csv"), row.names = F)


## mr_instrument
do.call(rbind.data.frame, lapply(mr_output.l, {
  function(l) l[["mr_het"]] }) ) %>%
  mutate(set=factor(Set, levels=sets, labels=names(sets)))  %>%
  write.csv(paste0(saveDir, "/", savePref,"_mr_instrument_", tag, ".csv"), row.names=F)

## mr_output.l 
saveRDS(mr_output.l, file=paste0(saveDir, "/", savePref, "_mr_output_", tag, ".rda"))

##EOF


