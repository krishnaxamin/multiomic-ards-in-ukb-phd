#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -N same_disease_ldsc
#$ -o ~/gwas/lava/logs/
#$ -t 1-68

# Run LDSC for all disease pairs

# LD panel (1KGP EUR) downloaded from: https://zenodo.org/records/8182036

mapfile -t diseases < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

# LDSC
module load miniforge
# activate ldsc_py3 venv and navigate to venv dir
conda activate ldsc_py3
cd ""

python ldsc.py \
--rg "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/munged/${disease}.sumstats.gz","~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/munged/${disease}.sumstats.gz" \
--ref-ld-chr ~/data/external/genomics/ldsc/eur_w_ld_chr/ \
--w-ld-chr ~/data/external/genomics/ldsc/eur_w_ld_chr/ \
--out "~/ch2_genomics/2.11 local_genetic_correlations/int_data/diseaseX_diseaseX_ldsc_results/${disease}-${disease}"
