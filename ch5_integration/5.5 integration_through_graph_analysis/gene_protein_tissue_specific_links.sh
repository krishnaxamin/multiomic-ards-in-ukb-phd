#!/bin/bash
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=64G
#SBATCH -J gene_protein_graph_analysis
#SBATCH -o gene_protein_graph_analysis/outputs/%x.o%A.%a
#SBATCH -a 1-39

module load miniforge
mamba activate phd_project_mamba

# use SLURM_ARRAY_TASK_ID for tissue specification within the py script

# gives custom error message if 'cd' fails
cd ~/cellular_network || { echo "cellular_network could not be found"; exit 1; }

# AD was a hotspot in the GWAS-proteomics pathway-based similarities
# python gene_protein_tissue_specific_links.py \
#   --gwas_disease "G30" \
#   --proteomics_disease "G30" \
#   --proteomics_algo "cox" \
#   --tissue_index $SLURM_ARRAY_TASK_ID

# exploring osteoarthritis being affected by metabolic dysfunction
# if necessary, multiple diseases must be provided in a single string and separated by '-', e.g. "M15-M16"
python gene_protein_tissue_specific_links.py \
  --gwas_disease "E11" \
  --proteomics_disease "M15-M16-M17-M18" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID

python gene_protein_tissue_specific_links.py \
  --gwas_disease "N18" \
  --proteomics_disease "M15-M16-M17-M18" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID

python gene_protein_tissue_specific_links.py \
  --gwas_disease "I25" \
  --proteomics_disease "M15-M16-M17-M18" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID

python gene_protein_tissue_specific_links.py \
  --gwas_disease "K76" \
  --proteomics_disease "M15-M16-M17-M18" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID

# exploring Parkinson's, osteoarthritis, cataracts, interstitial pulmonary disease link
# firth because proteomics-based relationships between these diseases only seen using prevalent diagnoses, not incident
python gene_protein_tissue_specific_links.py \
  --gwas_disease "M15-M16-M17-M18" \
  --proteomics_disease "G20" \
  --proteomics_algo "firth" \
  --tissue_index $SLURM_ARRAY_TASK_ID

python gene_protein_tissue_specific_links.py \
  --gwas_disease "M15-M16-M17-M18" \
  --proteomics_disease "H25-H26" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID
  
python gene_protein_tissue_specific_links.py \
  --gwas_disease "H25-H26" \
  --proteomics_disease "M15-M16-M17-M18" \
  --proteomics_algo "cox" \
  --tissue_index $SLURM_ARRAY_TASK_ID

python gene_protein_tissue_specific_links.py \
  --gwas_disease "G20" \
  --proteomics_disease "J84" \
  --proteomics_algo "firth" \
  --tissue_index $SLURM_ARRAY_TASK_ID
  
python gene_protein_tissue_specific_links.py \
  --gwas_disease "M15-M16-M17-M18" \
  --proteomics_disease "J84" \
  --proteomics_algo "firth" \
  --tissue_index $SLURM_ARRAY_TASK_ID