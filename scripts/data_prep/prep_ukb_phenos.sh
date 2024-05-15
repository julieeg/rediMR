#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=1:00:00
#$ -o ./reports

#$ -j y
#$ -cwd


source /broad/software/scripts/useuse
use R-4.1


Rscript ../scripts/data_prep/prep_ukb_phenos.R


#EOF

