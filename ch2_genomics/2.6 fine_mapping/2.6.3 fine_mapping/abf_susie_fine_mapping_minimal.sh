#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=100G
#$ -j y
#$ -N abf_susie_fine_mapping
#$ -o /data/WHRI-Phenogenomics/krishna/fine_mapping/saige/minimal/logs/
#$ -t 1-354

# Fine mapping of all disease-relevant regions not already covered in the colocalisation pipeline.
# Firstly, single-causal-variant ABF is run. 
# Then, multiple-causal-variant SuSiE is attempted. This will fail for some diseases due to GWAS signals being incompatible with the region's LD matrix.

module load R/4.4.1

Rscript ./abf_susie_fine_mapping.R ${SGE_TASK_ID} minimal
