#!/bin/bash


#$ -l h_vmem=5G
#$ -l h_rt=0:30:00
#$ -o reports/

#$ -cwd
#$ -j y



pheno=$1
tag=$2
ss_prefix=../data/processed/gwas/${pheno}_${tag}


source /broad/software/scripts/useuse
use R-4.1



# Merge chromosome-specific summary statistics
head -1 ${ss_prefix}_chr22.gwas > ${ss_prefix}.gwas

# choose TEST column
TEST=$(awk -v RS='\t' '/TEST/{print NR; exit}' ${ss_prefix}.gwas)

for i in {1..22}; do
	echo "Gathering summary status from ${ss_prefix}_chr${i}.gwas ..."
	awk -v TEST=$TEST '{ if ($TEST == "ADD") {print $0} }' OFS='\t' ${ss_prefix}_chr${i}.gwas >> ${ss_prefix}.gwas
done

head ${ss_prefix}.gwas

# Make gwas plots
Rscript ../scripts/gwas/postprocess_gwas.R ${pheno} ${tag}


## EOF
