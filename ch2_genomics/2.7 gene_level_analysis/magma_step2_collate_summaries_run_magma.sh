#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N magma_step2
#$ -o ~/gwas/magma/magma_step2_logs/
#$ -t 1-68

# MAGMA step 2: assemble all summary statistics, and run MAGMA

# set disease
mapfile -t disease_fields < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease_field=${disease_fields[${SGE_TASK_ID}]};

# set which GWAS adjustment to use
adjustment=$1

# assemble summary statistics
module load R;
# navigate to directory where the MAGMA executable is
# cd ~/gwas/magma
if [ -f "~/ch2_genomics/2.7 gene_level_analysis/gwas_summary_stats_for_magma/${disease_field}_magma_${adjustment}.tsv" ]; then
echo 'Summary statistics already collated.'
else 
Rscript magma_step2_collate_summaries.R ${disease_field} --adjustment ${adjustment}
fi

# run MAGMA
./magma \
--bfile ~/data/external/genomics/magma/g1000_eur synonyms=~/data/external/genomics/magma/g1000_eur.synonyms synonym-dup=drop-dup \
--gene-annot ~/data/external/genomics/magma/all_genotyped_imputed.genes.annot \
--pval "~/ch2_genomics/2.7 gene_level_analysis/gwas_summary_stats_for_magma/${disease_field}_magma_${adjustment}.tsv" use=snpid,pval duplicate=first  N=425223 \
--seed 42 \
--out ~/data/internal/genomics/magma/${adjustment}/gene_level_results/${disease_field}
