#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=24:00:00

#$ -j y
#$ -cwd


CHR=$SGE_TASK_ID

pheno=$1
covars="age sex gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10"


source /broad/software/scripts/useuse
use R-4.1
use UGER

reuse -q Anaconda3
source activate ../opt/bgen


scratch=/broad/hptmp/gervis

ukb_bgen_dir=/broad/ukbb/imputed_v3 #[ukb_imp_chr${chr}_v3.bgen
ukb_sample_dir=/humgen/florezlab/UKBB_app27892 #[ukb27892_imp_chrAUT_v3_s487395.sample]



# Genotype QC (maf>0.005 & info>0.4)
awk -v CHR=$CHR '{ if ( $6 > 0.005 && $8 > 0.4 ) { print $2 }  }' ${ukb_bgen_dir}/ukb_mfi_chr${CHR}_v3.txt > ../data/processed/gwas/snplist_chr${CHR}_maf005_imp04.txt

bgenix -g ${ukb_bgen_dir}/ukb_imp_chr${CHR}_v3.bgen -incl-rsids ../data/processed/gwas/snplist_chr${CHR}_maf005_imp04.txt > ${scratch}/chr${CHR}_sel_5k.bgen


# Format phenotype file using R
R --no-save <<EOF
library(tidyverse) ; library(data.table)
vars_to_select<-c("${pheno}", strsplit("${covars}", " ")[[1]])
fread("../data/processed/ukb_phenos_unrelated_EUR.csv", header=T) %>% mutate("FID"=id) %>% select(FID, IID=id, all_of(vars_to_select)) %>% write.table("../data/processed/gwas/${pheno}_chr${CHR}.txt", col.names=F, row.names=F, quote=F)
EOF

sed -i 1i"#FID IID ${pheno} ${covars}" ../data/processed/gwas/${pheno}_chr${CHR}.txt


echo "Done preparing phenotype file for analysis: ../data/processed/gwas/${pheno}_chr${CHR}.csv"
head ../data/processed/gwas/${pheno}_chr${CHR}.csv



## Run GWAS, by chromosome
../opt/plink2 \
	--bgen ${scratch}/chr${CHR}_sel_5k.bgen ref-first \
	--sample ${ukb_sample_dir}/ukb27892_imp_chrAUT_v3_s487395.sample \
	--id-delim \
	--glm \
	--pheno ../data/processed/gwas/${pheno}_chr${CHR}.txt \
	--pheno-name ${pheno} \
	--covar-name ${covars} \
	--covar-variance-standardize \
	--rm-dup force-first \
	--hwe 1e-12 \
	--geno 0.2 \
	--mind 0.2 \
	--memory 50000 \
	--threads 8 \
	--out ../data/processed/gwas/${pheno}_chr${CHR} \
&& rm ../data/processed/gwas/${pheno}_chr${CHR}.txt

## save logs
mv ../data/gwas/gwas_chr${CHR}_${pheno}_${covar_bases}.log ./logs/
rm ${scratch}/chr${CHR}_sel_5k.bgen

## END_OF_FILE
