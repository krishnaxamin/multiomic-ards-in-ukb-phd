#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N snp_qc_step7_hwe
#$ -o snp_qc_pipeline/step7_hwe/logs/
#$ -t 1-68

# calculate HWE stats with high missingness (>2%) variants and samples removed, and low MAF (<1%) variants removed. Step 7 of the SNP QC pipeline

# Use PLINK2 to identify which variants fail HWE < 10^-6 when using both female and male control populations for all chromosomes, as per this post https://groups.google.com/g/plink2-users/c/vA61ymOtZ5o/m/pBCKizX3AQAJ and this paper https://pubmed.ncbi.nlm.nih.gov/27071844/.
# This should output HWE stats for autosomes in one file and stats for cX in another file.

mapfile -t diseases < ../pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

module load plink;

plink2 \
--bfile /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--fam /data/WHRI-Phenogenomics/krishna/ukb_genotype_data/ukb_genotyping_col6_is_1.fam \
--keep control_eids/pan-ukbb-eur_${disease}-control_eids.tsv \
--remove snp_qc_pipeline/step4_high_missingness_samples_panukbb_eur.tsv \
--exclude snp_qc_pipeline/high_missingness_low_maf_variants_panukbb_eur.tsv \
--hardy midp \
--threads $NSLOTS \
--out snp_qc_pipeline/step7_hwe/${disease}_panukbb_eur