#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_models
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data5/UK_Biobank/16p12_general_population/3_second_hits
#SBATCH -o /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/13_go_models.log
#SBATCH -e /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/13_go_models.log

echo `date` starting job on $HOSTNAME

Rscript 13_go_logistic.R

echo `date` ending job
