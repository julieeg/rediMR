#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


## ReDiMR Inputs
pheno=$1
tag=$2
pheno_tag=${pheno}_${tag}

ssFile=../data/processed/gwas/${pheno_tag}.gwas
phenoFile=../data/processed/ukb_phenos_unrelated.rda


## Parameters
P=5e-8
pctBthhold=20
ANC=EUR

outDir=../data/processed/rediMR


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

awk -v p=$P '{ if($12 < p) {print $3} }' ${ssFile} > ${outDir}/${pheno_tag}_snps

# build bgen 
for i in {1..22}; do
	echo Builing bgen for chr ${i} ...
	bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${outDir}/${pheno_tag}_snps > ${scratch}/${pheno_tag}_chr${i}_snps.bgen
	done ; cat-bgen -g ${scratch}/${pheno_tag}_chr*_snps.bgen -og ${scratch}/${pheno_tag}_snps.bgen -clobber \
&& rm ${scratch}/${pheno_tag}_chr*_snps.bgen



echo Running LD clumping on ${pheno_tag}_snps and extracting dosage...

## run LD clumping in plink
../opt/plink2 \
--bgen ${scratch}/${pheno_tag}_snps.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--export A \
--freq \
--clump ${ssFile} \
--clump-p1 5e-8 \
--clump-r2 0.01 \
--clump-kb 250 \
--clump-field "P" \
--clump-snp-field "ID" \
--rm-dup force-first \
--memory 50000 \
--out ${outDir}/${pheno_tag}_snps \
&& mv ${outDir}/${pheno_tag}_snps.clumps ${outDir}/${pheno_tag}_clumps \
&& rm ${scratch}/${pheno_tag}_snps.bgen 


awk '{if(NR>1) {print $3} }' ${outDir}/${pheno_tag}_clumps > ${outDir}/${pheno_tag}_loci
n=$(wc -l < ${outDir}/${pheno_tag}_loci)



echo Gathering dosage data ...


## Extract dosage & afreq for loci, only
R --vanilla <<EOF 
library(tidyverse) ; library(data.table) 
loci<-c(fread("${outDir}/${pheno_tag}_loci", header=F))[[1]] 
fread("${outDir}/${pheno_tag}_snps.raw") %>% select(IID, contains(loci)) %>% fwrite("${outDir}/${pheno_tag}_loci.raw") 
fread("${outDir}/${pheno_tag}_snps.afreq") %>% filter(ID %in% loci) %>% fwrite("${outDir}/${pheno_tag}_loci.afreq") 
EOF



echo Done preparing genotype data for ${pheno}_${tag}: $n loci were identified. 



##########################################################
## Format summary stats & merge phenotypes/dosage files ##
##########################################################

echo Formatting summary stats and merging phenotype/dosage files in R.

R --vanilla <<EOF 
library(tidyverse) ; library(data.table) 

# Format summary stats file: writing ${datOutput}/${pheno_tag}_ssInput 
g<-fread("${outDir}/${pheno_tag}_loci.raw") 

fread("${ssFile}") %>% rename(CHR="#CHROM", SNP=ID, N=OBS_CT) %>% \
 mutate(EA=ifelse(A1==REF, REF, ALT), NEA=ifelse(A1==REF, ALT, REF), LOCI=ifelse(SNP %in% gsub("_[^_]*$", "", names(g)), 1, 0)) %>% \
 select(CHR, POS, SNP, EA, NEA, BETA, SE, P, N, LOCI) %>% \
 filter(P < 5e-8 & LOCI == 1) %>% fwrite("${outDir}/${pheno_tag}_ssInput.tmp.csv", row.names=F) 

# Merge phenotypes/dosage files: writing {outDir}/${pheno_tag}_datInput 
left_join(readRDS("${phenoFile}"), g %>% rename(id=IID) %>% mutate(id=as.character(id)), by = "id") %>% filter(ancestry == "${ANC}") %>% \
saveRDS("${outDir}/${pheno_tag}_datInput.tmp.rda") 
EOF



# Format datInput.rda file
Rscript ../scripts/rediMR/prep_datInput.R ${pheno} ${tag} 


############################
## Run ReDiMR program -v1 ##
############################

ssInput=${outDir}/${pheno_tag}_ssInput.csv 
datInput=${outDir}/${pheno_tag}_datInput.rda 
outDir=${outDir} 
pctBthhold=20 


Rscript ../scripts/rediMR/rediMR.R $pheno $tag $ssInput $datInput $pctBthhold $outDir



####################################
## Delete breadcrumbs from Step 1 ##
####################################

echo Deleting the following breadcrumb files from ${outDir}: \
${pheno_tag}_datInput.tmp.rda ${pheno_tag}_ssInput.tmp.csv ${pheno_tag}_loci ${pheno_tag}_loci.raw

rm ${outDir}/${pheno_tag}_*.raw

rm ${outDir}/${pheno_tag}_datInput.tmp.rda 
rm ${outDir}/${pheno_tag}_ssInput.tmp.csv 


echo --> You are now ready to proceed to ReDiMR Step 2-Mendelian Randomization. \
Good luck!


##EOF

