#!/bin/bash


#$ -l h_vmem=25G
#$ -l h_rt=2:00:00
#$ -o reports/

#$ -cwd
#$ -j y



testDir=../data/processed/testsnps
bgenDir=/broad/ukbb/imputed_v3
bgen_sample=/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample
scratch=/broad/hptmp/gervis


phenoFile=../data/processed/ukb_phenos_unrelated.rda


## Load resources
plink2=../../opt/plink2
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1



## Write file with 6 test snps
#touch ${testDir}/testsnps
#echo rs1421085 rs1726866 rs6690619 rs7619139 rs838133 rs9323534 > ${testDir}/testsnps



## Make bgen
for i in {1..22}; do
bgenix -g ${bgen_dir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${testDir}/testsnps > ${scratch}/testsnps_chr${i}.bgen
done ; cat-bgen -g ${scratch}/testsnps_chr*.bgen -og ${scratch}/testsnps.bgen -clobber \
&& rm ${scratch}/testsnps_chr*.bgen



# Export dosages
${plink2} \
--bgen ${scratch}/testsnps.bgen ref-first \
--sample $bgen_sample \
--export A \
--out ${testDir}/testsnps



#Combine dosage + phenotypes into ukb_analysis_test
R --vanilla <<EOF
# Merge phenotypes/dosage files: write {testDir}/${testsnps}_datInput 

EOF





