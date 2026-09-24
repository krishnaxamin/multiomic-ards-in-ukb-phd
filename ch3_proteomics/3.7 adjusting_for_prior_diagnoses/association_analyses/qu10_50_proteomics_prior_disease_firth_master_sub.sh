#!/bin/bash
#SBATCH -t 1:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=2G
#SBATCH -J qu10_50_proteomics_lifestyles_priors_firth_master_sub
#SBATCH -o /data/home/bty207/proteomics/association_analysis/firth_lifestyle_prior_disease/nested_subs_logs/%x.o%A.%a
#SBATCH -a 1-68

mapfile -t diseases < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SLURM_ARRAY_TASK_ID}]};

# test for presence of complete result - only submit if absent
if [ ! -f "~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/firth_results/${disease}_all_proteins.csv" ]; then
	sbatch ./qu10_50_proteomics_prior_disease_firth_nested_subs.sh --disease ${disease}
else
	echo Analysis for ${disease} already complete.
fi
