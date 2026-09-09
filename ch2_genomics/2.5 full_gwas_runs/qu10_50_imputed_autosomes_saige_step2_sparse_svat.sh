#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -j y
#$ -N qu10_50_imputed_autosomes_sparse_svat
#$ -o /data/WHRI-Phenogenomics/krishna/ukb_saige_sparse_svat/qu10_50_imputed_svat_logs/
#$ -t 1-1564

# Running SAIGE1.1.9 SVAT on the imputed autosomes data, for the 68 qu10_50 diseases. Submit is arrayed so that each job is one chromosome for one disease.

mapfile -t diseases_chrs < pan_ukbb_eur_qu10_50_chr_for_apocrita_arrayed_submit.txt;  # disease_chr file
disease_chr=${diseases_chrs[${SGE_TASK_ID}]};

disease="$(echo "$disease_chr" | cut -d'_' -f1)";
chr="$(echo "$disease_chr" | cut -d'_' -f2)";

mkdir "~/ch2_genomics/2.5 full_gwas_runs/step2_svat/${disease}";

# note that bgen files are treated alt-first here because these bgen files are made by PLINK, which outputs them as alt-first
apptainer run saige_1.1.9.sif step2_SPAtests.R \
--bgenFile=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen \
--bgenFileIndex=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen.bgi \
--sampleFile=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.sample \
--is_imputed_data=TRUE \
--minMAF=0.01 \
--minInfo=0.5 \
--chrom=${chr} \
--LOCO=FALSE \
--AlleleOrder=alt-first \
--SAIGEOutputFile="~/ch2_genomics/2.5 full_gwas_runs/step2_svat/${disease}/${disease}_c${chr}_pan-ukbb-eur.imputed.svat" \
--GMMATmodelFile="~/ch2_genomics/2.5 full_gwas_runs/step1_var_estimate/${disease}.rda" \
--varianceRatioFile="~/ch2_genomics/2.5 full_gwas_runs/step1_var_estimate/${disease}.varianceRatio.txt" \
--sparseGRMFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
--sparseGRMSampleIDFile="~/ch2_genomics/2.5 full_gwas_runs/step0_sparse_grm/pan_ukbb_eur_sparse_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
--is_output_moreDetails=TRUE \
--is_Firth_beta=TRUE \
--pCutoffforFirth=0.05 \
--is_fastTest=TRUE
