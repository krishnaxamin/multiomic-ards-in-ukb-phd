#!bin/bash
#$ -cwd
#$ -m bea
#$ -l h_rt=240:0:0
#$ -pe smp 16
#$ -l h_vmem=16G
#$ -l highmem
#$ -j y
#$ -N pan_ukbb_eur_saige_sparse_grm
#$ -o /data/scratch/bty207/saige_sparse_grm/saige_sparse_grm.log

# Shell script to get a Sparse GRM for the Pan-UKBB EUR cohort (Step 0)

module load plink;

# make bfile set with the Pan-UKBB EUR samples only. Write out the bfile set to the temp dir.
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--remove "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv" \
--output-chr 26 \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/autosomes-cX_for_sparse_grm;

# Step 0 - create sparse GRM
apptainer run saige_1.1.9.sif createSparseGRM.R \
--plinkFile=$TMPDIR/autosomes-cX_for_sparse_grm \
--minMAFforGRM=0.01 \
--maxMissingRateforGRM=0.02 \
--nThreads=$NSLOTS \
--outputPrefix="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm";
