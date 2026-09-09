#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N get_proteomics_qu10_50_prior_phenotypes
#$ -o get_proteomics_qu10_50_prior_phenotypes.log

# generate ARD proteomics prior phenotype files from the file downloaded off DNAnexus

module load miniforge
mamba activate phd_project_mamba

cd ~

python filter_dnanexus_phenotype_file.py \
--phenotype_dnanexus_file ~/data/external/genomics/array-genotyping_phenotype_file.csv \
--phenotype_covariate_file ~/data/internal/genomics/genomics_covars.csv \
--export_eids_filtered n \
--phenotype_set pan_ukbb_eur_qu10_50_disease_fields.txt \
-o ~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv
