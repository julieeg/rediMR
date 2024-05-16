#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


## ReDiMR Inputs
pheno=$1
tag=$2

ssFile=../data/processed/gwas/${pheno}.gwas
phenoFile=../data/processed/ukb_phenos_unrelated.rda


## Parameters
P=5e-8
pDelta=20
ANC=EUR

outDir=../data/processed/rediMR/${pheno}


## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1


genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis



#################################################
## Run LD clumping to extract significant loci ##
#################################################

mkdir -p ${outDir}


awk -v p=$P '{ if($12 < p) {print $3} }' ${ssFile} > ${outDir}/${pheno}.gwas.snps.${tag}

# build bgen 
for i in {1..22}; do
	echo Builing bgen for chr ${i} ...
	bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${outDir}/${pheno}.gwas.snps.${tag} > ${scratch}/${pheno}_chr${i}.gwas.snps.${tag}.bgen
	done ; cat-bgen -g ${scratch}/${pheno}_chr*.gwas.snps.${tag}.bgen -og ${scratch}/${pheno}.gwas.snps.${tag}.bgen -clobber \
&& rm ${scratch}/${pheno}_chr*.gwas.snps.${tag}.bgen



echo Running LD clumping on ${pheno}.gwas.snps.${tag} and extracting dosage...

## run LD clumping in plink
../opt/plink2 \
--bgen ${scratch}/${pheno}.gwas.snps.${tag}.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--export A \
--clump ${ssFile} \
--clump-p1 5e-8 \
--clump-r2 0.01 \
--clump-kb 250 \
--clump-field "P" \
--clump-snp-field "ID" \
--rm-dup force-first \
--memory 50000 \
--out ${outDir}/${pheno}.gwas.snps.${tag}

awk '{if(NR>1) {print $3} }' ${outDir}/${pheno}.gwas.snps.${tag}.clumps > ${outDir}/${pheno}.gwas.loci.${tag}
n=$(wc -l < ${outDir}/${pheno}.gwas.loci.${tag})



echo Gathering dosage data ...


## Extract dosage for loci, only
R --vanilla <<EOF 
library(tidyverse) ; library(data.table) 
loci<-c(fread("${outDir}/${pheno}.gwas.loci.${tag}", header=F))[[1]] 
fread("${outDir}/${pheno}.gwas.snps.${tag}.raw") %>% select(IID, contains(loci)) %>% fwrite("${outDir}/${pheno}.gwas.loci.${tag}.raw") 
EOF



echo Done preparing genotype data for $pheno: $n loci were identified. 



##########################################################
## Format summary stats & merge phenotypes/dosage files ##
##########################################################

echo Formatting summary stats and merging phenotype/dosage files in R.

R --vanilla <<EOF 
library(tidyverse) ; library(data.table) 

# Format summary stats file: writing ${datOutput}/${pheno}_ssInput 
g<-fread("${outDir}/${pheno}.gwas.loci.${tag}.raw") 

fread("${ssFile}") %>% rename(CHR="#CHROM", SNP=ID, N=OBS_CT) %>% \
 mutate(EA=ifelse(A1==REF, REF, ALT), NEA=ifelse(A1==REF, ALT, REF), LOCI=ifelse(SNP %in% gsub("_[^_]*$", "", names(g)), 1, 0)) %>% \
 select(CHR, POS, SNP, EA, NEA, BETA, SE, P, N, LOCI) %>% \
 filter(P < 5e-8 & LOCI == 1) %>% fwrite("${outDir}/${pheno}.${tag}_ssInput.tmp.csv", row.names=F) 

# Merge phenotypes/dosage files: writing {outDir}/${pheno}_datInput 
left_join(readRDS("${phenoFile}"), g %>% rename(id=IID) %>% mutate(id=as.character(id)), by = "id") %>% filter(ancestry == "${ANC}") %>% \
saveRDS("${outDir}/${pheno}.${tag}_datInput.tmp.rda") 
EOF



# Format datInput.rda file
Rscript ../scripts/rediMR/prep_datInput.R ${pheno} ${tag} \
&& rm ${outDir}/${pheno}.${tag}_datInput.tmp.rda



############################
## Run ReDiMR program -v1 ##
############################

ssInput=${outDir}/${pheno}_ssInput.csv 
datInput=${outDir}/${pheno}_datInput.rda 
outDir=${outDir} 
pctBchange=20 

tag=${tag}


Rscript ../scripts/rediMR/rediMR.R ${pheno} ${tag}



##EOF





