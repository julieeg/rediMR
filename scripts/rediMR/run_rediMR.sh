#!/bin/bash


#$ -l h_vmem=70G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


## ReDiMR Inputs
pheno=$1
tag=$2
covars=$3

pheno_tag=${pheno}_${tag}

redimrDir=../data/processed/rediMR/vCole/
ssInput=${redimrDir}/${pheno_tag}_ssInput.csv
datInput=${redimrDir}/${pheno_tag}_datInput.rda


# Default parameters
pctBthhold=20 
ANC=EUR


## Load resources
source /broad/software/scripts/useuse

use R-4.1



########################################
##  ~~~~~  Run ReDiMR program  ~~~~~  ##
########################################

Rscript ../scripts/rediMR/rediMR_vCole.R $pheno $tag $covars &&
echo You completed rediMR for $pheno, adjusting for $covarSet covariates.



##EOF

