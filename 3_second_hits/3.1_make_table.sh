#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=ukb_table
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/UK_Biobank/16p12_general_population/3_second_hits
#SBATCH -o /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/3_make_table/%a.log
#SBATCH -e /data5/UK_Biobank/16p12_general_population/3_second_hits/logs/3_make_table/%a.log
#SBATCH --array 1-24

echo `date` starting job on $HOSTNAME

chrom=$(head -n $SLURM_ARRAY_TASK_ID analysis_files/chromosomes.list | tail -n 1)

echo $chrom

python 3_secondhit_table.py $chrom

echo `date` ending job
