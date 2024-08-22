#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=5:00:00
#$ -o reports/

#$ -j y
#$ -cwd



pheno=$1
tag=$2


pheno_tag=${pheno}_${tag}
ANC=EUR

bgen_dir=/broad/ukbb/imputed_v3
bgen_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
scratch=/broad/hptmp/gervis

prs_dir=../data/processed/prs
redimr_dir=../data/processed/rediMR




## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1




#################################################
## Build bgen file for gwas loci ##
#################################################

## make array 
declare -a snpsets=("all" "ref" "unref")

for SET in "${snpsets[@]}" ; do
	
	# Build bgen file for each loci set
	echo Building bgen file for the ${SET} SNP set of ${pheno_tag} ...
	for i in {1..22}; do
	bgenix -g ${bgen_dir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${redimr_dir}/${pheno_tag}_loci > ${scratch}/${pheno_tag}_chr${i}_loci.bgen
	done ; cat-bgen -g ${scratch}/${pheno_tag}_chr*_loci.bgen -og ${scratch}/${pheno_tag}_loci.bgen -clobber \
	&& rm ${scratch}/${pheno_tag}_chr*_loci.bgen ;
done


for SET in "${snpsets[@]}" ; do

	# Build prs
	echo Running PRS for the ${SET} SNP set of ${pheno_tag} ...
	../opt/plink2 \
	--bgen ${scratch}/${pheno_tag}_loci.bgen ref-first \
	--sample $bgen_sample \
	--maf 0.005 \
	--score ${prs_dir}/${pheno_tag}_${SET}_prsInput \
	--rm-dup force-first \
	--out ${prs_dir}/${pheno_tag}_${SET}

done





plink --bgen