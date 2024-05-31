#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=16:00:00
#$ -o reports/

#$ -j y
#$ -cwd



ANC=$1


ukb_bgen_dir=/broad/ukbb/imputed_v3 #[ukb_imp_chr${chr}_v3.bgen
ukb_sample_dir=/humgen/florezlab/UKBB_app27892 #[ukb27892_imp_chrAUT_v3_s487395.sample]

scratch=/broad/hptmp/gervis


source /broad/software/scripts/useuse
use R-4.1

reuse -q Anaconda3
source activate ../opt/bgen




## Prepare phenotype file (format as #FID IID; as tsv)
R --no-save <<EOF
library(dplyr) ; library(data.table) 
fread("../data/processed/ukb_phenos_unrelated.csv", header=T) %>% filter(ancestry == "${ANC}") %>% \
mutate(FID=id, IID=id) %>% \
relocate(FID, IID, .before=id) %>% \
mutate(across(colnames(.), function(x) gsub(" ", "_", x))) %>% rename('#FID'=FID) %>% \
write.table("../data/processed/gwas/ukb_phenos_unrelated_${ANC}_gwas.txt", col.names=T, row.names=F, quote=F)
EOF



## Prepare genotype data

# filter SNPs by maf>0.005 & info>0.4
awk -v CHR=$CHR '{ if ( $6 > 0.005 && $8 > 0.4 ) { print $2 }  }' ${ukb_bgen_dir}/ukb_mfi_chr${CHR}_v3.txt > ../data/processed/gwas/snplist_chr${CHR}_maf005_imp04.txt

# Make bgen files by chrom
bgenix -g ${ukb_bgen_dir}/ukb_imp_chr${CHR}_v3.bgen -incl-rsids ../data/processed/gwas/snplist_chr${CHR}_maf005_imp04.txt > ${scratch}/chr${CHR}_sel.bgen



#EOF 

