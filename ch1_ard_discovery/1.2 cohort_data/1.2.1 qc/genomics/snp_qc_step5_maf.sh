#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -m bea
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N snp_qc_step5_maf
#$ -o snp_qc_pipeline/step5_maf/snp_qc_step5_maf.log

# calculate minor allele frequency with high missingness (>2%) variants and samples removed. Step 5 of the SNP QC pipeline

module load plink

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--not-chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--freq counts \
--threads $NSLOTS \
--out snp_qc_pipeline/step5_maf/autosomes_panukbb_eur;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_male_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--freq counts \
--threads $NSLOTS \
--out snp_qc_pipeline/step5_maf/cX_panukbb_eur_male;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--chr X \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep pan-ukbb-eur_female_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--exclude snp_qc_pipeline/step2_high_missingness_variants_panukbb_eur.tsv \
--freq counts \
--threads $NSLOTS \
--out snp_qc_pipeline/step5_maf/cX_panukbb_eur_female;