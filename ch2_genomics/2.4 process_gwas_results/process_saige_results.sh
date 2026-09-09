#!bin/bash
#$ -cwd
#$ -m bea
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=50G
#$ -j y
#$ -N process_saige_results
#$ -o process_saige_results.log

module load miniforge
mamba activate phd_project_mamba

cd ~/gwas

python process_saige_results.py  --adjustment "$1"