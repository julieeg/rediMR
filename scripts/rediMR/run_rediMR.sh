#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


## ReDiMR Inputs
pheno=$1

ssFile=../data/processed/gwas/${pheno}.gwas
phenoFile=../data/processed/ukb_phenos_unrelated.rda

P=5e-8
pDelta=0.2 #$2

outDir=../data/processed/rediMR/${pheno}


## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis



#################################################
## Run LD clumping to extract significant loci ##
#################################################

mkdir -p ${outDir}

awk -v p=$P '{ if($12 < p) {print $3} }' ${ssFile} > ${outDir}/${pheno}.gwas.snps.${P}

# build bgen 
for i in {1..22}; do
	echo Builing bgen for chr ${i} ...
	bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${outDir}/${pheno}.gwas.snps.${P} > ${scratch}/${pheno}_chr${i}.gwas.snps.${P}.bgen
	done ; cat-bgen -g ${scratch}/${pheno}_chr*.gwas.snps.${P}.bgen -og ${scratch}/${pheno}.gwas.snps.${P}.bgen -clobber \
&& rm ${scratch}/${pheno}_chr*.gwas.snps.${P}.bgen


echo Converting ${pheno}.gwas.snps.${P}.bgen to ${pheno}.gwas.snps.${P}.bfile ...

# convert bgen to bfile (plink2) 
../opt/plink2 --bgen ${scratch}/${pheno}.gwas.snps.${P}.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--make-bed \
--memory 50000 \
--rm-dup force-first \
--out ${scratch}/${pheno}.gwas.snps.${P}


echo Running LD clumping on ${pheno}.gwas.snps.${P} ...

## run LD clumping in plink
../opt/plink \
--bfile ${scratch}/${pheno}.gwas.snps.${P} \
--clump ${ssFile} \
--clump-p1 ${P} \
--clump-r2 0.01 \
--clump-kb 250 \
--clump-field "P" \
--clump-snp-field "ID" \
--out ${outDir}/${pheno}.gwas.loci.${P}

awk '{if(NR>1) {print $3} }' ${outDir}/${pheno}.gwas.loci.${P}.clumped > ${outDir}/${pheno}.gwas.loci.${P}
n=$(wc -l < ${outDir}/${pheno}.gwas.loci.${P})


echo Extracting dosage for $n loci for $pheno


## Extract dosage for loci in plink2 
../opt/plink2 --bgen ${scratch}/${pheno}.gwas.snps.${P}.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--extract ${outDir}/${pheno}.gwas.loci.${P} \
--export A \
--rm-dup force-first \
--memory 50000 \
--out ${outDir}/${pheno}.gwas.loci.${P} \
&& rm ${scratch}/${pheno}.gwas.snps.${P}.* \
&& rm ${outDir}/${pheno}.gwas.snps.${P} \
&& rm ${outDir}/${pheno}.gwas.loci.${P}.nosex


echo Done gathering genotype data.

##########################################################
## Format summary stats & merge phenotypes/dosage files ##
##########################################################

echo Formatting summary stats and merging phenotype/dosage files in R.

source /broad/software/scripts/useuse
use R-4.1


R --vanilla <<EOF 
library(tidyverse) ; library(data.table) 

# Format summary stats file: writing ${datOutput}/${pheno}_ssInput 
g<-fread("${outDir}/${pheno}.gwas.loci.${P}.raw") 
fread("${ssFile}") %>% rename(CHR="#CHROM", SNP=ID, N=OBS_CT) %>% \
 mutate(EA=ifelse(A1==REF, REF, ALT), NEA=ifelse(A1==REF, ALT, REF), LOCI=ifelse(SNP %in% gsub("_[^_]*$", "", names(g)), 1, 0)) %>% \
 select(CHR, POS, SNP, EA, NEA, BETA, SE, P, N, LOCI) %>% \
 filter(P < 5e-8) %>% fwrite("${outDir}/${pheno}_ssInput.csv", row.names=F) 

# Merge phenotypes/dosage files: writing {outDir}/${pheno}_datInput 
left_join(readRDS("${phenoFile}"), g %>% rename(id=IID), by = "id") %>% \
saveRDS("${outDir}/${pheno}_datInput_tmp.rda") 
EOF


# Format datInput file ; save as rda
Rscript ../scripts/rediMR/prep_datInput.R ${pheno} \
&& rm ${outDir}/${pheno}_datInput.tmp


############################
## Run ReDiMR program -v1 ##
############################

ssInput=${outDir}/${pheno}_ssInput.csv 
datInput=${outDir}/${pheno}_datInput.rda 
outDir=${outDir} 
pctBchange=20 


Rscript ../scripts/rediMR/rediMR.R ${pheno}



##EOF





