#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 240:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH -J gene_protein_in_diff_clusters
#SBATCH -o gene_protein_distances/clustering/outputs/%x.o%A.%a
#SBATCH -a 1-3

module load miniforge
mamba activate phd_project_mamba

# use SLURM_ARRAY_TASK_ID to do arrays of the three possible edge options. also use this to set the logname/jobname in the py script
# run each script twice for each of Firth and Cox proteomics links

edges=('interacts_with' 'coexpresses_with' 'both')

cd ~/cellular_network

python gene_protein_in_diff_clusters.py --gene_edges_to_use ${edges[$((SLURM_ARRAY_TASK_ID-1))]} --proteomics_algo_to_use 'firth'
python gene_protein_in_diff_clusters.py --gene_edges_to_use ${edges[$((SLURM_ARRAY_TASK_ID-1))]} --proteomics_algo_to_use 'cox'
