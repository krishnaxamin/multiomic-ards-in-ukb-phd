#!/bin/bash
#SBATCH -t 1:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=10G
#SBATCH -J proteomics_prior_disease_before_ards_pc_effects
#SBATCH -o prior_disease_info/effect_on_phenotype_logs/%x.o%A.%a
#SBATCH -a 1-68

module load R/4.4.1;

echo 'Before-ARDs PCs effect on prior (prevalent) phenotype.'
Rscript prior_disease_pc_before_ards_effects_on_phenotype.R ${SLURM_ARRAY_TASK_ID}
