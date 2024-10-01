#!/bin/bash
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -cwd
#$ -j y


exposure=$1 #oilyfish_QT
outcome=$2 #tg
sumstats=$3 #std_mrdat_oilyfish_GCST90239664_TG_Graham_GRCh37.csv
covarset=$4 #confounders1
tag=$5 #vCole

sumstatsfile=../data/sumstats/standardized_shared_from_KSJ/${sumstats}


saveDir=../data/processed/rediMR/vCole/${exposure}_${outcome}
mkdir -p $saveDir


## Load resources
source /broad/software/scripts/useuse
use R-4.1


## Run rediMR step2 
Rscript --no-save ../scripts/rediMR/rediMR_step2.R ${exposure} ${outcome} ${sumstatsfile} ${covarset} ${tag} ${saveDir} 

##EOF


