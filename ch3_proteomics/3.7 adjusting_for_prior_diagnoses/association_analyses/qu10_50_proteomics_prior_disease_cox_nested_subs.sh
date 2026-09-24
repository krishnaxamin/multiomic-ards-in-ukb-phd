#!/bin/bash
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=16G
#SBATCH -J qu10_50_proteomics_lifestyles_priors_cox
#SBATCH -o /data/home/bty207/proteomics/association_analysis/cox_lifestyle_prior_disease/nested_subs_logs/%x.o%A.%a
#SBATCH -a 1-2920

# association analysis
module load R/4.4.1;
Rscript ./qu10_50_proteomics_prior_disease_cox.R ${SLURM_ARRAY_TASK_ID} "$@";

# script to check if all disease-proteins runs have completed for a given disease
# if a disease-protein run is hanging (i.e., hasn't completed after 24h), this step will still think the disease task is incomplete, so collation must be manual in that case
Rscript ./qu10_50_proteomics_prior_disease_cox_collate_nested_subs.R ${SLURM_ARRAY_TASK_ID} "$@";

echo 'Association complete.'
