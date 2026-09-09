#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -pe smp 1
#$ -l h_vmem=100G
#$ -j y
#$ -N ld_calcs_for_assoc_blocks
#$ -o /data/WHRI-Phenogenomics/krishna/ld_calcs/ld_mats/logs/
#$ -t 1-331

# Calculate r values for LD between variants in predefined regions. Part of the pipelines performing fine mapping, colocalisation, and regional analysis/plotting of associated variants.

# load in bgen module (bgen v1.1.7)
module load miniforge
mamba activate bgen_env

# load in predefined regions and isolate which region we're looking at in this task
mapfile -t regions < blocks_for_ld_calcs.csv;
IFS=',', read -r chr start end <<< "${regions[${SGE_TASK_ID}]}"

# filter imputed data from the correct chromosome and with the pan-ukbb-eur samples to the range of the predefined region
bgenix -g ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.bgen -incl-range ${chr}:${start}-${end} > /data/scratch/bty207/ld_calcs/task${SGE_TASK_ID}_imputed_data.bgen

# create index file for new bgen file
bgenix -g ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data.bgen -index

# unload BGEN; load PLINK2.0 for reading BGEN input
mamba deactivate
module load plink

# convert the filtered imputed data to a bed/bim/fam dataset so it can be merged with the genotype data. also make snpIDs to be of form chr:pos:ref:alt (set-all-var-ids)
plink2 \
--bgen ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data.bgen ref-last \
--sample ~/data/external/genomics/ukb_imputed_data/pan_ukbb_eur/c${chr}_pan-ukbb-eur_info-filtered.sample \
--set-all-var-ids @:#:\$r:\$a \
--new-id-max-allele-len 1000 \
--make-bed \
--threads $NSLOTS \
--out ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data

# filter genotype data to correct samples and correct chromosome and range, and make snpIDs to be of form chr:pos:ref:alt
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--remove "~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv" \
--chr ${chr} \
--from-bp ${start} \
--to-bp ${end} \
--set-all-var-ids @:#:\$r:\$a \
--new-id-max-allele-len 1000 \
--threads $NSLOTS \
--make-bed \
--out ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_genotyped_data

# unload PLINK2.0, load PLINK1.9
module unload plink;
module load plink/1.9-beta6.27-gcc-12.2.0;

# test for if there is genotyped data for this region. If so read it in and merge (--bmerge) in the binary PLINK fileset for the filtered imputed data.
# if not, then read in the binary PLINK fileset for the filtered imputed data only.
# calculate r values, with output being a square matrix (--r square).
# also write out of the list of snpIDs used to calculate r values. 
# --keep-allele-order is necessary to prevent PLINK1.9 from changing the allele order so that the minor allele is the effect allele
#  see this issue in susieR: https://github.com/stephenslab/susieR/issues/148
if [ -e "~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_genotyped_data.bed" ]; then
plink \
--bfile ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_genotyped_data \
--bmerge ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data \
--write-snplist \
--keep-allele-order \
--r square \
--threads $NSLOTS \
--out ~/data/external/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}
else 
plink \
--bfile ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data \
--write-snplist \
--keep-allele-order \
--r square \
--threads $NSLOTS \
--out ~/data/external/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}
fi

if [ -e "~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_genotyped_data.bed" ]; then
rm ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_genotyped_data*
fi;
rm ~/data/internal/genomics/ld_calcs/task${SGE_TASK_ID}_imputed_data*
if [ -e "~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.bed" ]; then 
rm ~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.bed
fi;
if [ -e "~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.bim" ]; then
rm ~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.bim
fi;
if [ -e "~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.fam" ]; then
rm ~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_${chr}_${start}_${end}.fam
fi;
