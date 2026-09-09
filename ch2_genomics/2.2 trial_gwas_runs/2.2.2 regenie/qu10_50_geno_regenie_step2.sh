#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -j y
#$ -N qu10_50_geno_regenie_step2
#$ -o /data/scratch/bty207/regenie/step2/qu10_50_geno_logs/
#$ -t 1-69

# Running REGENIE3.3 Step 2 on the genotyped data, for the 69 qu10_50 diseases.

module load plink;

mapfile -t diseases < ~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

# REGENIE does not recognise 25 as a valid chromosome code - replace with XY/PAR1/PAR2. 'MT' option below will also replace 23 with 'X'.
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--exclude "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/variants_to_exclude/${disease}_autosomes-cX-vars_to_exclude.tsv" \
--not-chr X \
--output-chr MT \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/${disease}-for-step2

module unload plink;

module load miniconda;

mamba activate regenie_env;

# mkdir /data/scratch/bty207/regenie/step2/${disease};

regenie \
--step 2 \
--bed $TMPDIR/${disease}-for-step2 \
--phenoFile ~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/${disease}_pheno-file_REGENIE.tsv \
--phenoCol ${disease} \
--covarFile ~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/${disease}_covars_REGENIE.tsv \
--covarColList p21003_i0,p31,p54_i0,p22000,PC{1:20} \
--catCovarList p31,p54_i0,p22000 \
--maxCatLevels 110 \
--pred "~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step1/${disease}_pan-ukbb-eur.step1_pred.list" \
--bt \
--bsize 1000 \
--threads $NSLOTS \
--firth \
--approx \
--firth-se \
--pThresh 0.05 \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step2/${disease}/${disease}_pan-ukbb-eur.geno.step2"
