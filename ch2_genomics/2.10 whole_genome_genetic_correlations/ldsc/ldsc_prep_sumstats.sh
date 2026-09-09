#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -N ldsc_prep_sumstats
#$ -o ~/gwas/ldsc/logs/
#$ -t 1-68

# Script to prep summary stats into LDSC-compatible format

# alleles file for --merge-alleles in munge_sumstats.py: https://zenodo.org/records/7773502 (HapMap3 SNP list)

mapfile -t diseases < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};
echo 'prepping for ${disease}'

# own python script to collate and prep sumstats for munging
module load miniforge
conda activate phd_project_mamba

cd "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc"
python prep_sumstats_for_ldsc.py --disease ${disease}
echo "prep complete for $disease"

# munging
conda deactivate
# activate ldsc_py3 venv and navigate to venv dir
conda activate ldsc_py3
cd ""

python ./munge_sumstats.py \
--sumstats "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/collated/${disease}.tsv" \
--frq MAF \
--merge-alleles /data/home/bty207/gwas/ldsc/w_hm3.snplist \
--out "~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/munged/${disease}"
echo "munge complete for $disease"

