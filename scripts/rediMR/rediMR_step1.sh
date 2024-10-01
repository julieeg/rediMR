#!/bin/bash
#$ -l h_vmem=30G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


exposure=$1 #oilyfish_QT
outcome=$2 #tg
sumstats=$3 #std_mrdat_oilyfish_GCST90239664_TG_Graham_GRCh37.csv
covarset=$4 #confounders1
tag=$5 #vCole

phenofile=../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits_09292024.txt 
covars_gwas="age sex gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10"

saveDir=../data/processed/rediMR/vCole/${exposure}_${outcome}
mkdir -p $saveDir

## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1


genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis


##############################################
## Gather GWAS SNPs from Summary Stats File ##
##############################################

sumstatsfile=../data/sumstats/standardized_shared_from_KSJ/${sumstats}

## First if sumstats is a csv, convert to txt --> Then, grab list of SNPs
R --vanilla <<EOF
library(data.table) ; library(tidyverse)
ss=fread("${sumstatsfile}")
if(endsWith("${sumstatsfile}", ".csv")) {
  "sumstats is a csv --> converting to txt"
  ss %>% fwrite(gsub("csv", "txt", "${sumstatsfile}"), sep=" ") } else{
    "sumstats is a txt --| No conversion required."
  } ; ss %>% select("snp.exposure") %>% fwrite("${saveDir}/${exposure}_${outcome}.snpsInput", sep=" ", col.names=F)
EOF

snp_list=${saveDir}/${exposure}_${outcome}.snpsInput


############################################
## Prepare individual-level genotype data ##
############################################

# build bgen file with input snps (.snpsInput)
for i in {1..22}; do
echo Builing bgen for chr ${i} ...
bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${snp_list} > ${scratch}/${exposure}_${outcome}_snps_chr${i}_${tag}.bgen
done ; cat-bgen -g ${scratch}/${exposure}_${outcome}_snps_chr*_${tag}.bgen -og ${scratch}/${exposure}_${outcome}_snps_${tag}.bgen -clobber \
&& rm ${scratch}/${pheno}_${outcome}_snps_chr*_${tag}.bgen


## Get dosage & AF
../opt/plink2 \
--bgen ${scratch}/${exposure}_${outcome}_snps_${tag}.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--extract ${snp_list} \
--export A \
--freq \
--rm-dup force-first \
--memory 50000 \
--out ${saveDir}/${exposure}_${outcome}_snpsInput \
&& rm ${scratch}/${exposure}_${outcome}_snps_${tag}.bgen 


## Run ReDiMR - Step 1: SNP Refinement using Covariate Adjustment
Rscript --no-save ../scripts/rediMR/rediMR_step1.R $exposure $outcome $sumstatsfile $covarset $tag $saveDir

