#!/bin/bash
#$ -l h_vmem=60G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


pheno=$1 #oilyfish_QT
tag=$2 #vCole
r2=$3
#pheno_tag=${pheno}_${tag}
covars=confounders
ssFile=../data/processed/gwas/${pheno}.gwas #oilyfish_QT_vCole.gwas
phenoFile=../data/processed/ukb_phenos_unrelated_EUR_withJC_diet_traits.txt #$3
redimrDir=../data/processed/rediMR/vCole #$4


## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1


genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis
gwasDir=../data/processed/gwas


############################################
## Prepare individual-level genotype data ##
############################################

## Run LD clumping & Get SNP dosage ----------------

cat ${ssFile} | awk '{ if($10 < 5e-8) {print $1} }' > ${ssFile}.snps

# build bgen 
for i in {1..22}; do
echo Builing bgen for chr ${i} ...
bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${ssFile}.snps > ${scratch}/${pheno_tag}_chr${i}_snps.bgen
done ; cat-bgen -g ${scratch}/${pheno_tag}_chr*_snps.bgen -og ${scratch}/${pheno_tag}_snps.bgen -clobber \
&& rm ${scratch}/${pheno_tag}_chr*_snps.bgen


echo Running LD clumping on ${pheno_tag}_snps and extracting dosage...


## run LD clumping in plink
../opt/plink2 \
--bgen ${scratch}/${pheno_tag}_snps.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--clump ${ssFile} \
--clump-p1 5e-8 \
--clump-r2 ${r2} \
--clump-kb 500 \
--clump-field "P" \
--clump-snp-field "ID" \
--rm-dup force-first \
--memory 50000 \
--out ${ssFile}_${r2}_.snps

awk '{if(NR>1) {print $3} }' ${ssFile}_${r2}.snps.clumps > ${ssFile}_${r2}.loci
n=$(wc -l < ${ssFile}_${r2}.loci)


## get dosage & allele frequency for loci
../opt/plink2 \
--bgen ${scratch}/${pheno_tag}_snps.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--extract ${ssFile}_${r2}.loci \
--export A \
--freq \
--rm-dup force-first \
--memory 50000 \
--out ${ssFile}_${r2}.loci


#&& mv ${ssFile}.clumps ${ssFile}/${pheno_tag}_clumps #\
rm ${scratch}/${pheno_tag}_snps.bgen 


## Format & merge individual-level data in R ----------------
Rscript --no-save ../scripts/rediMR/prep_rediMR_vCole.R $pheno $tag $covars $r2


#EOF



