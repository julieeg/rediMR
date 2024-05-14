#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=0:30:00
#$ -o reports/

#$ -cwd
#$ -j y



pheno=$1

ss_path=../data/processed/gwas/${pheno}


ukb_bgen_dir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis


reuse -q Anaconda3
source activate ../opt/bgen



# Merge chromosome-specific summary statistics
head -1 ${ss_path}/${pheno}_chr${CHR}.gwas > ${ss_path}/${pheno}_merged.gwas
for CHR in {i..22}; do
	echo "${ss_path}/${pheno}_chr${CHR}.gwas ..."
	tail -n +2 ${ss_path}/${pheno}_chr${CHR}.gwas >> ${ss_path}/${pheno}_merged.gwas
done



# Subset to significant loci (P<0.08)
awk '{ if ($7 == "ADD" && $12 < 5e-8) { print $3 } }' ${ss_path}/${pheno}_merged.gwas > ${ss_path}/${pheno}.loci


# Create plinkfile for ld_clumping
for CHR in {1..22}; do
bgenix -g /broad/ukbb/imputed_v3/ukb_imp_chr${CHR}_v3.bgen -incl-rsids ${ss_path}/${pheno}.loci > ${scratch}/${pheno}_chr${CHR}.loci.bgen
done

cat-bgen -g ${scratch}/${pheno}_chr${CHR}loci.bgen -og ${scratch}/${pheno}_loci.bgen -clobber



# Run ld clumping
../opt/plink2 --bgen ${scratch}/$${pheno}_loci.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--make-bed \
--memory 5000 \
--rm-dup force-first \
--out ${scratch}/${pheno}.loci


../opt/plink \
--bfile ${scratch}/${pheno}.loci \
--clump ../data/processed/gwas/${pheno}_merged.gwas \
--clump-p1 5e-8 \
--clump-r2 0.002 \
--clump-kb 250 \
--clump-field "P" \
--clump-snp-field "ID" \
--out ../data/processed/gwas/${pheno}.loci




# Format summary stats in R for downstream analysis
Rscript ../scripts/postprocess_gwas.R ${pheno}



## EOF



