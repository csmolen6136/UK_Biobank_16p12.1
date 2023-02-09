#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_models
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data5/UK_Biobank/16p12_general_population/1_icd10_codes
#SBATCH -o /data5/UK_Biobank/16p12_general_population/1_icd10_codes/logs/5_models.log
#SBATCH -e /data5/UK_Biobank/16p12_general_population/1_icd10_codes/logs/5_models.log

echo `date` starting job on $HOSTNAME

Rscript 5_logistic_model.R

echo `date` ending job
