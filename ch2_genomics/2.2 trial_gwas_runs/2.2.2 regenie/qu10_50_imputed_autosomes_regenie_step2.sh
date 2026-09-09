#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -j y
#$ -N qu10_50_imputed_regenie_step2
#$ -o /data/scratch/bty207/regenie/step2/qu10_50_imputed_logs/
#$ -t 1-1587

# Running REGENIE3.3 on the imputed autosomes data, for the 69 qu10_50 diseases. Submit is arrayed so that each job is one chromosome for one disease.

mapfile -t diseases_chrs < "~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/qu10_50_chr_for_apocrita_arrayed_submit.txt";  # disease_chr file
disease_chr=${diseases_chrs[${SGE_TASK_ID}]};

disease="$(echo "$disease_chr" | cut -d'_' -f1)";
chr="$(echo "$disease_chr" | cut -d'_' -f2)";

module load miniconda;

mamba activate regenie_env;

regenie \
--step 2 \
--bgen ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen \
--sample ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.sample \
--ref-first \
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
--chr ${chr} \
--minINFO 0.5 \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step2/${disease}/${disease}_c${chr}_pan-ukbb-eur.imputed.step2"
