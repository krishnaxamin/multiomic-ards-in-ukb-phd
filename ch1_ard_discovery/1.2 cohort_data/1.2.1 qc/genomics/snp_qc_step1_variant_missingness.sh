#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -j y
#$ -N snp_qc_step1_variant_missingness
#$ -o results/snp_qc_pipeline/step1_variant_missingness/snp_qc_step1_variant_missingness.log

# calculate variant missingness. uses the pan-ukbb EUR population in the genomic cohort. Step 1 of the SNP QC pipeline

module load plink

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--not-chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_eids.tsv \
--missing variant-only \
--threads $NSLOTS \
--out results/snp_qc_pipeline/step1_variant_missingness/autosomes_panukbb_eur;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_male_eids.tsv \
--missing variant-only \
--threads $NSLOTS \
--out results/snp_qc_pipeline/step1_variant_missingness/cX_panukbb_eur_male;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_female_eids.tsv \
--missing variant-only \
--threads $NSLOTS \
--out results/snp_qc_pipeline/step1_variant_missingness/cX_panukbb_eur_female;