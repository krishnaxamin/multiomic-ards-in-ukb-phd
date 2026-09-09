#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N snp_qc_step3_sample_missingness
#$ -o snp_qc_pipeline/step3_sample_missingness/snp_qc_step3_sample_missingness.log

# calculate sample missingness with high missingness (>2%) variants removed. Step 3 of the SNP QC pipeline

module load plink

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--not-chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_eids.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--missing sample-only \
--threads $NSLOTS \
--out snp_qc_pipeline/step3_sample_missingness/autosomes_panukbb_eur;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_male_eids.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--missing sample-only \
--threads $NSLOTS \
--out snp_qc_pipeline/step3_sample_missingness/cX_panukbb_eur_male;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_female_eids.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--missing sample-only \
--threads $NSLOTS \
--out snp_qc_pipeline/step3_sample_missingness/cX_panukbb_eur_female;