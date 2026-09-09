#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -j y
#$ -N qu10_50_geno_sparse_svat
#$ -o /data/scratch/bty207/saige_sparse_svat/qu10_50_geno_svat_logs/
#$ -t 1-68

# Perform SAIGE1.1.9 SVAT on the PLINK data for the 68 qu10_50 diseases.

module load plink;

mapfile -t diseases < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

# make bfile set with the variants and samples to keep. Write out the bfile set to the temp dir.
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--remove "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv" \
--exclude "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/variants_to_exclude/${disease}_autosomes-cX-vars_to_exclude.tsv" \
--output-chr MT \
--make-bed \
--threads $NSLOTS \
--out "~/ch2_genomics/2.5 full_gwas_runs/svat_temp_dir/${disease}-for-svat"

# chr_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X XY);

# loop through the c1-22
# for i in ${chr_list[@]};
# make bfiles just for that chromosome. Note that the chr coding here defaults to 'MT' (numeric for autosomes, then 'X','Y','XY'/'PAR1'/'PAR2' as necessary)
# do plink2 \
# --bfile /data/scratch/bty207/svat_temp_dir/${disease}-for-svat \
# --chr ${i} \
# --make-bed \
# --threads $NSLOTS \
# --out /data/scratch/bty207/svat_temp_dir/c${i}_${disease}-for-svat;

mkdir "~/ch2_genomics/2.5 full_gwas_runs/step2_svat/${disease}";

# run SAIGE SVAT
apptainer run saige_1.1.9.sif step2_SPAtests.R \
--bedFile="~/ch2_genomics/2.5 full_gwas_runs/svat_temp_dir/${disease}-for-svat.bed" \
--bimFile="~/ch2_genomics/2.5 full_gwas_runs/svat_temp_dir/${disease}-for-svat.bim" \
--famFile="~/ch2_genomics/2.5 full_gwas_runs/svat_temp_dir/${disease}-for-svat.fam" \
--LOCO=FALSE \
--AlleleOrder=alt-first \
--SAIGEOutputFile="~/ch2_genomics/2.5 full_gwas_runs/step2_svat/${disease}/${disease}_pan-ukbb-eur.geno.svat" \
--GMMATmodelFile="~/ch2_genomics/2.5 full_gwas_runs/step1_var_estimate/${disease}.rda" \
--varianceRatioFile="~/ch2_genomics/2.5 full_gwas_runs/step1_var_estimate/${disease}.varianceRatio.txt" \
--sparseGRMFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
--sparseGRMSampleIDFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
--is_output_moreDetails=TRUE \
--is_Firth_beta=TRUE \
--pCutoffforFirth=0.05 \
--is_fastTest=TRUE

# done;
rm "~/ch2_genomics/2.5 full_gwas_runs/svat_temp_dir/${disease}*"