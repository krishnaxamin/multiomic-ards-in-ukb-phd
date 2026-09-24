#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=32G
#SBATCH -J olink_cossims_vs_all
#SBATCH -o gene_protein_distances/embeddings/outputs/%x.o%A

module load miniforge
mamba activate phd_project_mamba

# gives custom error message if 'cd' fails
cd ~/cellular_network || { echo "cellular_network could not be found"; exit 1; }

python olink_cossims_vs_all.py
