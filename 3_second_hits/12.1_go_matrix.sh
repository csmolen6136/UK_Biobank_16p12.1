#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_models
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data5/UK_Biobank/16p12_general_population/3_second_hits
#SBATCH -o /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/12_go_matrix/%a.log
#SBATCH -e /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/12_go_matrix/%a.log
#SBATCH --array 1-100
#SBATCH --exclude sarah,qingyu

echo `date` starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME

python 12_count_tables.py $SLURM_ARRAY_TASK_ID

echo `date` ending job
