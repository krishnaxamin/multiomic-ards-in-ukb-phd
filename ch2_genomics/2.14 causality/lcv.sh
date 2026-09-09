#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -j y
#$ -N lcv
#$ -o ./logs/
#$ -t 1-2278

mapfile -t disease_pairs < ./disease_pairs.txt;
disease_pair=${disease_pairs[${SGE_TASK_ID}]};

module load R;

Rscript lcv.R ${disease_pair}