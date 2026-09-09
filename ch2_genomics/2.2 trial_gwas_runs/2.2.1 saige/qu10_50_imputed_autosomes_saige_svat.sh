#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -j y
#$ -N qu10_50_imputed_autosomes_svat
#$ -o /data/scratch/bty207/saige_svat/qu10_50_imputed_svat_logs/
#$ -t 1-1587

# Running SAIGE1.1.9 SVAT on the imputed autosomes data, for the 69 qu10_50 diseases. Submit is arrayed so that each job is one chromosome for one disease.

mapfile -t diseases_chrs < "~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/qu10_50_chr_for_apocrita_arrayed_submit.txt";  # disease_chr file
disease_chr=${diseases_chrs[${SGE_TASK_ID}]};

disease="$(echo "$disease_chr" | cut -d'_' -f1)";
chr="$(echo "$disease_chr" | cut -d'_' -f2)";


apptainer run saige_1.1.9.sif step2_SPAtests.R \
--bgenFile=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen \
--bgenFileIndex=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen.bgi \
--sampleFile=~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.sample \
--is_imputed_data=TRUE \
--minMAF=0.01 \
--minInfo=0.5 \
--chrom=${chr} \
--LOCO=True \
--AlleleOrder=ref-first \
--SAIGEOutputFile="~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/saige_svat/${disease}/${disease}_imputed_c${chr}_pan-ukbb-eur.imputed.svat" \
--GMMATmodelFile="~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/saige_grm/${disease}.rda" \
--varianceRatioFile="~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/saige_grm/${disease}.varianceRatio.txt" \
--is_output_moreDetails=TRUE \
--is_Firth_beta=TRUE \
--pCutoffforFirth=0.05 \
