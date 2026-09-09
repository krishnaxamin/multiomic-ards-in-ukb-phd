#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=60G
#$ -j y
#$ -N qu10_50_regenie_step1
#$ -o /data/WHRI-Phenogenomics/krishna/ukb_regenie_step1/qu10_50_logs/
#$ -t 1-69

# Shell script to run REGENIE3.3 Step 1 on the qu10_50 diseases, on Apocrita.

module load plink;

mapfile -t diseases < ~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

# REGENIE does not recognise 25 as a valid chromosome code - replace with XY/PAR1/PAR2. 'MT' option below will also replace 23 with 'X'.
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--extract ~/data/internal/genomics/ld_pruning/${disease}/${disease}_ld_pruning.prune.in \
--output-chr MT \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/autosomes-cX_${disease}-for-grm;

module unload plink;

module load miniconda;

mamba activate regenie_env;

regenie \
--step 1 \
--bed $TMPDIR/autosomes-cX_${disease}-for-grm \
--phenoFile ~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/${disease}_pheno-file_REGENIE.tsv \
--phenoCol ${disease} \
--covarFile ~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/${disease}_covars_REGENIE.tsv \
--covarColList p21003_i0,p31,p54_i0,p22000,PC{1:20} \
--catCovarList p31,p54_i0,p22000 \
--maxCatLevels 110 \
--bt \
--bsize 1000 \
--threads $NSLOTS \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step1/${disease}_pan-ukbb-eur.step1"
