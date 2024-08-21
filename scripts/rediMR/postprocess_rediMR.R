# Postprocess rediMR


#############################################
## ~~~~~~ Set up & input parameters ~~~~~~ ##
#############################################

# load required packages
lapply(c("tidyverse", "data.table", "parallel", "paletteer", "RColorBrewer",
         "ggpubr", "R3port", "tinytex"),  
       library, character.only = TRUE)


# load pantry file with stored parameters, basic functions & ggplot templates
source("../scripts/pantry.R")


# ==========================
## Load command arguments
# ==========================

# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1] #diet phenotype
tag <- args[2] #version "tag"
covarGroup <- args[3]
baseCovars= <- args[4]

redimrDir <- "../data/processed/rediMR/"

sets <- covarSetGroups$geneticPCs$Sets


#############################################
## ~~~~~~ Load & Format rediMR data ~~~~~~ ##
#############################################

do.call(rbind.data.frame, lapply(sets, function(set) {
  file <- paste0(pheno, "_", tag, "_", set, "_", baseCovars)
  list(all <- read.csv(paste0(redimrDir, "/", file, "_tabBchangeAll.csv")))[[1]]
  })) %>%
  mutate(gPCs=as.numeric(gsub("Top", "", gsub("geneticPCs", "", Covar)))) %>%
  mutate(snp=gsub("snp","",gsub("[.]",":",snp)),
         Covar=factor(Covar, levels=covarSetGroups$geneticPCs$Labels)) %>%
  fwrite(paste0("../data/processed/rediMR/", pheno, "_", tag, "_", covarGroup, "_", baseCovars, "_tabBchangeAll.csv"))






