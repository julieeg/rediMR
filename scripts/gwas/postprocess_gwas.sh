#!/bin/bash


#$ -l h_vmem=5G
#$ -l h_rt=0:30:00
#$ -o reports/

#$ -cwd
#$ -j y


pheno=$1
ss_prefix=../data/processed/gwas/${pheno}


source /broad/software/scripts/useuse
use R-4.1



# Merge chromosome-specific summary statistics
head -1 ${ss_prefix}_chr22.gwas > ${ss_prefix}.gwas
for i in {1..22}; do
	echo "Gathering summary status from ${ss_prefix}_chr${i}.gwas ..."
	awk '{ if ($7 == "ADD") {print $0} }' ${ss_prefix}_chr${i}.gwas >> ${ss_prefix}.gwas
done



# Make gwas plots
Rscript ../scripts/gwas/postprocess_gwas.R ${pheno}



## EOF
