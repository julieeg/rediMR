#$ -l h_vmem=30G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y



pheno=$1
tag=$2

pheno_tag=${pheno}_${tag}

ssFile=../data/processed/gwas/${pheno_tag}.gwas
phenoFile=../data/processed/ukb_phenos_unrelated.rda
redimrDir=../data/processed/rediMR


## Default parameters
ANC=EUR
gwasP=5e-8


## Load resources
source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate ../../opt/bgen

use R-4.1


genoDir=/broad/ukbb/imputed_v3 
scratch=/broad/hptmp/gervis



#################################################
## Run LD clumping to extract significant loci ##
#################################################

awk -v p=$gwasP '{ if($12 < p) {print $3} }' ${ssFile} > ${redimrDir}/${pheno_tag}_snps

# build bgen 
for i in {1..22}; do
echo Builing bgen for chr ${i} ...
bgenix -g ${genoDir}/ukb_imp_chr${i}_v3.bgen -incl-rsids ${redimrDir}/${pheno_tag}_snps > ${scratch}/${pheno_tag}_chr${i}_snps.bgen
done ; cat-bgen -g ${scratch}/${pheno_tag}_chr*_snps.bgen -og ${scratch}/${pheno_tag}_snps.bgen -clobber \
&& rm ${scratch}/${pheno_tag}_chr*_snps.bgen



echo Running LD clumping on ${pheno_tag}_snps and extracting dosage...

## run LD clumping in plink
../opt/plink2 \
--bgen ${scratch}/${pheno_tag}_snps.bgen ref-first \
--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
--export A \
--freq \
--clump ${ssFile} \
--clump-p1 5e-8 \
--clump-r2 0.01 \
--clump-kb 250 \
--clump-field "P" \
--clump-snp-field "ID" \
--rm-dup force-first \
--memory 50000 \
--out ${redimrDir}/${pheno_tag}_snps \
&& mv ${redimrDir}/${pheno_tag}_snps.clumps ${redimrDir}/${pheno_tag}_clumps \
&& rm ${scratch}/${pheno_tag}_snps.bgen 


awk '{if(NR>1) {print $3} }' ${redimrDir}/${pheno_tag}_clumps > ${redimrDir}/${pheno_tag}_loci
n=$(wc -l < ${redimrDir}/${pheno_tag}_loci)

echo LD pruning completed: $n independent SNPs identified at r2 of 0.01 with 250 kb clumping window




##END_OF_SYNTAX

