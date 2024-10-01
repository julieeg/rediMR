#prep_rediMR_dat.R

library(tidyverse) ; library(data.table) 


#command args
args <- commandArgs(trailingOnly = T)
pheno=args[1] #"oilyfish_QT"
tag=args[2] #"vCole"
covars=args[3] #"confounders"
r2=args[4] #LD clumping r2

pheno_tag=paste0(pheno,"_",tag)

ssFile=paste0("../data/processed/gwas/",pheno_tag,".gwas_", r2) #$2
phenoFile="../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits.txt" #$3
redimrDir="../data/processed/rediMR/vCole" #$4


source("../scripts/pantry.R")


# =========================================================
# Prepare individual-level 
# ==========================================================
cat("Loading individual-level genotype files")

## Load list of loci
loci_raw <- fread(paste0(ssFile, ".loci.raw")) %>% 
  rename_with(., ~gsub(":",".", .)) %>%
  rename_with(., ~ifelse(startsWith(., "rs") | . == "IID", ., paste0("id",.)))
  
# Prepare summary stats file
loci <- fread(paste0(ssFile, ".loci"), header=F)[[1]]
loci <- ifelse(startsWith(loci, "rs"), loci, paste0("id", gsub(":",".",loci)))

loci_ss <- fread(ssFile) %>% rename(SNP=ID) %>% 
  mutate_at("SNP", ~ifelse(startsWith(., "rs"), ., paste0("id", gsub(":",".",.)))) %>%
  filter(P < 5e-8) %>% 
  filter(SNP %in% loci) %>%
  select(CHR, POS=BP, SNP, ALT, REF, EAF, BETA, SE, P, INFO) %>%
  mutate(ID=paste0(CHR,":",POS,"_",ALT,"_",REF)) 

# Prepare loci.afreq file 
#fread(paste0(ssFile, ".snps.afreq")) %>% filter(ID %in% loci) %>% fwrite(paste0(ssFile, ".loci.afreq"))



# ====================================================================
# Prepare confounders: Add numeric versions of descriptive factors
# ====================================================================

dat <- fread(phenoFile) 
 
## Get confounders
covars_id <- dat %>% select(
  id, ancestry, age, sex, ac, paste0("gPC", 1:20),
  all_of(covarSets$confounders$Covars), all_of(covarSets$confounders_num$Covars),
  ldl, hdl, tg, chol, glu, crp, water)

## Get diet phenotypes (from JC)
phenos_id <- dat %>% select(id, oilyfish_QT, alch_glasspermonth_QT, bread_type_BIN)


# =========================================================
# Merge individual-level data into rediMR datInput
# ==========================================================

## Merge phenotype / genotype data ----------- 
left_join(phenos_id, covars_id, by="id") %>% 
  left_join(loci_raw %>% rename(id=IID) %>% mutate(id=as.integer(id)), by="id") %>% 
  filter(ancestry == "EUR") %>%
  saveRDS(paste0(redimrDir, "/", pheno_tag, "_datInput.rda")) 

loci_ss %>% fwrite(paste0(redimrDir, "/", pheno_tag, "_ssInput.csv"), row.names=F)

##EOF


