#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -o reports/

#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID

pheno=$1
covars="age sex gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10"



ukb_sample_dir=/humgen/florezlab/UKBB_app27892
scratch=/broad/hptmp/gervis


source /broad/software/scripts/useuse
use R-4.1

reuse -q Anaconda3
source activate ../opt/bgen



# Format phenotype file in R
R --no-save <<EOF
library(tidyverse) ; library(data.table)
vars_to_select<-c("${pheno}", strsplit("${covars}", " ")[[1]])
fread("../data/processed/ukb_phenos_unrelated_EUR.csv", header=T) %>% mutate("FID"=id) %>% select(FID, IID=id, all_of(vars_to_select)) %>% write.table("../data/processed/gwas/${pheno}_chr${CHR}.txt", col.names=F, row.names=F, quote=F)
EOF

sed -i 1i"#FID IID ${pheno} ${covars}" ../data/processed/gwas/${pheno}_chr${CHR}.txt


echo "Done preparing phenotype file for analysis: ../data/processed/gwas/${pheno}_chr${CHR}.txt"
head ../data/processed/gwas/${pheno}_chr${CHR}.txt



## Run GWAS, by chromosome
../opt/plink2 \
	--bgen ${scratch}/chr${CHR}_sel_test.bgen ref-first \
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
	--out ../data/processed/gwas/chr${CHR}_test \
&& mv ../data/processed/gwas/chr${CHR}_test.${pheno}.glm.linear ../data/processed/gwas/${pheno}_chr${CHR}_test.gwas



#EOF




