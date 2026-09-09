#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -N ldsc
#$ -o ~/gwas/ldsc/logs/
#$ -t 1-2278

# Run LDSC for all disease pairs

# LD panel (1KGP EUR) downloaded from: https://zenodo.org/records/8182036

mapfile -t disease_pairs < "../disease_pairs.txt";
disease_pair=${disease_pairs[${SGE_TASK_ID}]};

IFS='-' read -r disease1 disease2 <<< "$disease_pair"

# LDSC
module load miniforge
# activate ldsc_py3 venv and navigate to venv dir
conda activate ldsc_py3
cd ""

python ldsc.py \
--rg "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/munged/${disease1}.sumstats.gz","~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/munged/${disease2}.sumstats.gz" \
--ref-ld-chr ~/data/external/genomics/ldsc/eur_w_ld_chr/ \
--w-ld-chr ~/data/external/genomics/ldsc/eur_w_ld_chr/ \
--out "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/results/${disease_pair}"
