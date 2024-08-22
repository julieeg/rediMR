#$ -l h_vmem=30G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y



pheno=$1
tag=$2

pheno_tag=${pheno}_${tag}

ssFile=../data/processed/gwas/${pheno_tag}.gwas
phenoFile=../data/processed/ukb_phenos_unrelated.rda
redimrDir=../data/processed/rediMR


## Default parameters
ANC=EUR
gwasP=5e-8


## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1


genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis



#################################################################
## Prepare dosage; Format summary stats; Merge pheno/geno data ##
#################################################################

# use R

R --vanilla <<EOF 

library(tidyverse) ; library(data.table) 


## 1. Gathering dosage & allele frequency data ...

loci<-c(fread("${redimrDir}/${pheno_tag}_loci", header=F))[[1]] 
fread("${redimrDir}/${pheno_tag}_snps.raw") %>% select(IID, contains(loci)) %>% fwrite("${redimrDir}/${pheno_tag}_loci.raw") 
fread("${redimrDir}/${pheno_tag}_snps.afreq") %>% filter(ID %in% loci) %>% fwrite("${redimrDir}/${pheno_tag}_loci.afreq") 

# --> Done preparing genotype data for ${pheno}_${tag}: $n loci were identified.



# 2. Format summary stats file: writing ${datOutput}/${pheno_tag}_ssInput 

g<-fread("${redimrDir}/${pheno_tag}_loci.raw") 
fread("${ssFile}") %>% rename(CHR="#CHROM", SNP=ID, N=OBS_CT) %>% \
mutate(EA=ifelse(A1==REF, REF, ALT), NEA=ifelse(A1==REF, ALT, REF), LOCI=ifelse(SNP %in% gsub("_[^_]*$", "", names(g)), 1, 0)) %>% \
select(CHR, POS, SNP, EA, NEA, BETA, SE, P, N, LOCI) %>% \
filter(P < 5e-8 & LOCI == 1) %>% fwrite("${redimrDir}/${pheno_tag}_ssInput.tmp.csv", row.names=F) 



# 3. Merging phenotypes/dosage files: writing {redimrDir}/${pheno_tag}_datInput 

left_join(readRDS("${phenoFile}"), g %>% rename(id=IID) %>% mutate(id=as.character(id)), by = "id") %>% filter(ancestry == "${ANC}") %>% \
saveRDS("${redimrDir}/${pheno_tag}_datInput.tmp.rda") 

# DONE! The following files are ready: ${pheno_tag}_datInput.tmp.rda & ${pheno_tag}_ssInput.tmp.csv")


EOF


## Run prep_datInput.R
Rscript ../scripts/rediMR/prep_rediMR.R ${pheno} ${tag}

echo Done preparing ${pheno_tag}_datInput.rda and ${pheno_tag}_ssInput.csv files for rediMR!


## Delete breadcrumb files
rm $redimrDir/${pheno_tag}_*Input.tmp*
rm $redimrDir/${pheno_tag}*raw
rm $redimrDir/${pheno_tag}*snps.*


#EOF

