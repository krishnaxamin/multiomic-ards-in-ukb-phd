#!/bin/bash
#SBATCH -t 240:0:0
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --mem-per-cpu=50G
#SBATCH -J get_prior_disease_pcs_proteomics
#SBATCH -o /data/home/bty207/proteomics/prior_disease_info/get_prior_disease_pcs_proteomics.log

# get principal components representing ~50% variance explained by prior common-unisex disease information in the proteomics cohort
# prior disease information = diagnoses blood draw

cd "~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses"

module load miniforge
mamba activate phd_project_mamba

echo 'Getting PCs from diagnoses before blood draw and before ARD development'

python ./get_prior_disease_pcs.py \
--phenotype_dnanexus_file ~/data/internal/proteomics/proteomics_pan_ukbb_eur_phenotypes.csv \
--phenotype_covariate_file ~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv \
--background_phenotype_set ~/data/internal/phenotype_coding/common_unisex_disease_date_fields.csv \
--diagnosis_time_relation_to_ards before \
--diagnosis_time_relation_to_blood_draw before \
--cumulative_variance_path "~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/plots/common-unisex_disease_PCAs_cum_var_before_blood_before_ards" \
--pca_export_dir ~/data/internal/proteomics/prior_disease_info/principal_components_before_blood_before_ards \

echo 'Getting PCs from diagnoses before ARD development, agnostic to blood draw'

python ./get_prior_disease_pcs.py \
--phenotype_dnanexus_file ~/data/internal/proteomics/proteomics_pan_ukbb_eur_phenotypes.csv \
--phenotype_covariate_file ~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv \
--background_phenotype_set ~/data/internal/phenotype_coding/common_unisex_disease_date_fields.csv \
--diagnosis_time_relation_to_ards before \
--cumulative_variance_path "~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/plots/common-unisex_disease_PCAs_cum_var_before_ards" \
--pca_export_dir ~/data/internal/proteomics/prior_disease_info/principal_components_before_ards \
