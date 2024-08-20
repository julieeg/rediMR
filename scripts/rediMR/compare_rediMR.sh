#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y



pheno=$1
tag=$2
group=$3

pheno_tag=${pheno}_${tag}


## Load resources
source /broad/software/scripts/useuse

use R-4.1



## Run compare_rediMR.R

Rscript ../scripts/rediMR/compare_rediMR.R $pheno $tag $group


#EOF

