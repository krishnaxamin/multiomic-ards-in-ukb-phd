#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N get_metabolomics_qu10_50_prior_phenotypes
#$ -o get_metabolomics_qu10_50_prior_phenotypes.log

# generate ARD metabolomics prior phenotype files from the file downloaded off DNAnexus

module load miniforge
mamba activate phd_project_mamba

cd ~

python filter_dnanexus_phenotype_file.py \
--phenotype_dnanexus_file ~/data/external/metabolomics/metabolomics_phenotype_file_off_dnanexus.csv \
--phenotype_covariate_file ~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_covars.csv \
--export_eids_filtered n \
--phenotype_set pan_ukbb_eur_qu10_50_disease_fields.txt \
--phenotype_time prior \
-o ~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv
