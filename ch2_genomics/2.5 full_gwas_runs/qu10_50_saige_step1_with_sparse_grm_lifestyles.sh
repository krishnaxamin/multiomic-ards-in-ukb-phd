#!bin/bash
#$ -cwd
#$ -m bea
#$ -l h_rt=24:0:0
#$ -pe smp 1
#$ -l h_vmem=60G
#$ -j y
#$ -N qu10_50_saige_var_estimate_lifestyles
#$ -o /data/scratch/bty207/saige_var_estimate_lifestyles/qu10_50_logs/
#$ -t 1-68

# Shell script to estimate the variance ratio using a sparse GRM (Step 1) on the 68 qu10_50 diseases, on Apocrita.

module load plink;

mapfile -t diseases < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease=${diseases[${SGE_TASK_ID}]};

# make bfile set with the variants and samples to keep. Write out the bfile set to the temp dir. The previous SAIGE GRM run did not use the --remove option as no samples to remove due to missingness.
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--remove "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv" \
--extract ~/data/internal/genomics/ld_pruning/${disease}/${disease}_ld_pruning.prune.in \
--output-chr 26 \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/autosomes-cX_${disease}-for-grm;

apptainer run saige_1.1.9.sif step1_fitNULLGLMM.R \
--sparseGRMFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
--sparseGRMSampleIDFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
--useSparseGRMtoFitNULL=TRUE \
--skipVarianceRatioEstimation=FALSE \
--plinkFile=$TMPDIR/autosomes-cX_${disease}-for-grm \
--phenoFile=~/data/internal/genomics/pheno_covar_files_for_gwas/SAIGE/${disease}_pheno-file_SAIGE_lifestyles.tsv \
--phenoCol=${disease} \
--traitType=binary \
--covarColList=p21003_i0,p31,p54_i0,p23099_i0,p189,p884_i0,p1160_i0,p1289_i0,p1299_i0,p1309_i0,p1319_i0,p1438_i0,p1458_i0,p1528_i0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--qCovarColList=p31,p54_i0 \
--sampleIDColinphenoFile=IID \
--minMAFforGRM=0 \
--maxMissingRateforGRM=1 \
--nThreads=$NSLOTS \
--includeNonautoMarkersforVarRatio=TRUE \
--outputPrefix="~/ch2_genomics/2.5 full_gwas_runs/step1_var_estimate_lifestyles/${disease}" \
--IsOverwriteVarianceRatioFile=TRUE;
