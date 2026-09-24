#!/bin/bash
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=16G
#SBATCH -J gene_protein_embedding_distances
#SBATCH -o gene_protein_distances/embeddings/outputs/%x.o%A.%a
#SBATCH -a 1-234

module load miniforge
mamba activate phd_project_mamba

# use SLURM_ARRAY_TASK_ID for tissue specification within the py script

# gives custom error message if 'cd' fails
cd ~/cellular_network || { echo "cellular_network could not be found"; exit 1; }

python gene_protein_embedding_distances.py --task_index $SLURM_ARRAY_TASK_ID
