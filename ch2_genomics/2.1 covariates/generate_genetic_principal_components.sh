#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 3
#$ -l h_vmem=20G
#$ -j y
#$ -N ld_pca_arrayed
#$ -o /data/scratch/bty207/ld_pca_arrayed_logs/
#$ -t 1-68

# Combine LD pruning and PCA into one script for arrayed submit to Apocrita.

module load plink/1.9-170906;

mapfile -t disease_fields < ../pan_ukbb_eur_qu10_50_disease_fields.txt;
disease_field=${disease_fields[${SGE_TASK_ID}]};

mkdir /data/scratch/bty207/ld_pruning/${disease_field};

# LD pruning
plink \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--exclude variants_to_exclude/${disease_field}_autosomes-cX-vars_to_exclude.tsv \
--keep pan-ukbb-eur_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--indep-pairwise 500 1 0.2 \
--ld-xchr 1 \
--threads $NSLOTS \
--out /data/scratch/bty207/ld_pruning/${disease_field}/${disease_field}_ld_pruning;

module unload plink/1.9-170906;

module load plink;

mkdir /data/scratch/bty207/pca/${disease_field};

# PCA
plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--not-chr X \
--extract /data/scratch/bty207/ld_pruning/${disease_field}/${disease_field}_ld_pruning.prune.in \
--keep pan-ukbb-eur_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--pca allele-wts 20 approx vzs \
--threads $NSLOTS \
--out /data/scratch/bty207/pca/${disease_field}/${disease_field}_pca