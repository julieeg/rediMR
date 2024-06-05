#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o ./reports

#$ -j y
#$ -cwd


pheno=$1
tag=$2



source /broad/software/scripts/useuse
use R-4.1


Rscript ../scripts/mr/twosamplMR.R $pheno $tag


#EOF
