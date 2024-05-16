#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -o reports/

#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID


pheno=$1
phenoFile=../data/processed/ukb_phenos_unrelated_EUR.csv
covars="age sex gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10"

tag=$2


ukb_sample_dir=/humgen/florezlab/UKBB_app27892
scratch=/broad/hptmp/gervis



## Run GWAS, by chromosome
../opt/plink2 \
	--bgen ${scratch}/chr${CHR}_sel.bgen ref-first \
	--sample ${ukb_sample_dir}/ukb27892_imp_chrAUT_v3_s487395.sample \
	--id-delim \
	--glm \
	--pheno ${phenoFile} \
	--pheno-name ${pheno} \
	--covar-name ${covars} \
	--covar-variance-standardize \
	--rm-dup force-first \
	--hwe 1e-12 \
	--geno 0.2 \
	--mind 0.2 \
	--memory 50000 \
	--threads 8 \
	--out ${scratch}/${pheno}_chr${CHR}.${tag} \
&& mv ${scratch}/${pheno}_chr${CHR}.${tag}.${pheno}.glm.linear ../data/processed/gwas/${pheno}_chr${CHR}.${tag}.gwas \
&& mv ${scratch}/${pheno}_chr${CHR}.${tag}.log ../data/processed/gwas/${pheno}_chr${CHR}.${tag}.log 




#EOF


