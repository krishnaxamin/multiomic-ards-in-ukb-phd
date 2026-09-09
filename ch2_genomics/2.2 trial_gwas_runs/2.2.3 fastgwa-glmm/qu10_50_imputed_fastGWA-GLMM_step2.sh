#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -m bea
#$ -pe smp 2
#$ -l h_vmem=15G
#$ -j y
#$ -N qu10_50_imputed_fastGWA_step2
#$ -o /data/scratch/bty207/gcta/step2/testing/qu10_50_imputed_logs/

# testing for GCTA1.94.1 (fastGWA-GLMM) step 2, using c19 imputed data.

mapfile -t diseases_chrs < "~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/qu10_50_chr_for_apocrita_arrayed_submit.txt";  # disease_chr file
disease_chr=${diseases_chrs[${SGE_TASK_ID}]};

disease="$(echo "$disease_chr" | cut -d'_' -f1)";
chr="$(echo "$disease_chr" | cut -d'_' -f2)";

mkdir "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step2/${disease}"

./gcta-1.94.1 \
--bgen ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/${chr}_pan-ukbb-eur_info-filtered.bgen \
--sample ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/${chr}_pan-ukbb-eur_info-filtered.sample \
--grm-sparse "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}_grm_sparse" \
--fastGWA-mlm-binary \
--pheno ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/${disease}_pheno-file_GCTA.tsv \
--qcovar ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/${disease}_quant_covars_GCTA.tsv \
--covar ~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/categorical_covars_GCTA.tsv \
--threads $NSLOTS \
--maf 0.01 \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step2/${disease}/${disease}_c${chr}.imputed.assoc"
