#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -j y
#$ -N qu10_50_geno_fastGWA_step2
#$ -o /data/scratch/bty207/gcta/step2/qu10_50_geno_logs/
#$ -t 1-69

# GCTA1.94.1 (fastGWA-GLMM) step 2 for qu10_50 diseases.

mapfile -t diseases < ~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_disease_fields.txt;;
disease=${diseases[${SGE_TASK_ID}]};

module load plink;

plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--exclude "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/variants_to_exclude/${disease}_autosomes-cX-vars_to_exclude.tsv" \
--not-chr X \
--output-chr MT \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/${disease}_autosomes-cX_pan_ukbb_eur;

mkdir "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step2/${disease}";

./gcta-1.94.1 \
--bfile $TMPDIR/${disease}_autosomes-cX_pan_ukbb_eur \
--grm-sparse "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}_grm_sparse" \
--fastGWA-mlm-binary \
--pheno ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/${disease}_pheno-file_GCTA.tsv \
--qcovar ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/${disease}_quant_covars_GCTA.tsv \
--covar ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/categorical_covars_GCTA.tsv \
--threads $NSLOTS \
--maf 0 \
--geno 1 \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step2/${disease}/${disease}.geno.assoc"
