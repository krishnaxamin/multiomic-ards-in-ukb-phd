#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=60G
#$ -j y
#$ -N abf_susie_coloc
#$ -o ./logs/
#$ -t 1-95

# For each region where there are common associated variants between two diseaes, perform ABF (single causal variant) fine mapping (for each disease) and colocalisation. 
# Then, try SuSiE (multiple causal variant) fine mapping (for each disease), creating SuSiE objects for each disease. This will fail for some diseases due to GWAS signals being incompatible with the region's LD matrix. 
# If SuSiE fine mapping is successful for both disease, perform colocalisation using the SuSiE objects. 
# This is the final part of the pipeline colocalising genetic signals for ARDs, after 
#  (1) finding the regions containing common association variants 
#  (2)calculating the LD matrices for those regions. 

module load R/4.4.1

Rscript ~/gwas/colocalisations/saige/abf_susie_coloc.R ${SGE_TASK_ID}
