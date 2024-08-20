#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


## ReDiMR Inputs
pheno=$1
tag=$2
covarSet=$3
outcome=$4

pheno_tag=${pheno}_${tag}

redimrDir=../data/processed/rediMR
mrDir=../data/processed/MR

ssInput=${redimrDir}/${pheno_tag}_ssInput.csv 
datInput=${redimrDir}/${pheno_tag}_datInput.rda 


# Default parameters
pctBthhold=20 
ANC=EUR


## Load resources
source /broad/software/scripts/useuse

use R-4.1



######################################
##  ~~~~~  Run Two-Sample MR ~~~~~  ##
######################################

Rscript ../scripts/MR/MR.R $pheno $tag $covarSet $outcome $ssInput $datInput $pctBthhold $mrDir &&
echo Done with two sample MR!
  
  
#EOF

