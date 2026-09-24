#!/bin/bash
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=2G
#SBATCH -J qu10_50_proteomics_lifestyles_priors_cox_master_sub
#SBATCH -o /data/home/bty207/proteomics/association_analysis/cox_lifestyle_prior_disease/master_sub.log

# /data/home/bty207/proteomics/association_analysis/cox_lifestyle_prior_disease/nested_subs_logs/%x.o%A.%a

mapfile -t diseases < ../pan_ukbb_eur_qu10_50_disease_fields.txt;
#disease=${diseases[${SLURM_ARRAY_TASK_ID}]};

# test for presence of complete result - only submit if absent
for disease in ${diseases[@]}; do
  if [ ! -f "~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/cox_results/${disease}_all_proteins.csv" ]; then
    sbatch ./qu10_50_proteomics_prior_disease_cox_nested_subs.sh --disease ${disease}
    sleep 5m
  else
    echo Analysis for ${disease} already complete.
  fi
done
