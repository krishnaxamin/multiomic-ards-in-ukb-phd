#!bin/bash
#$ -cwd
#$ -l h_rt=240:0:0
#$ -m bea
#$ -pe smp 3
#$ -l h_vmem=10G
#$ -j y
#$ -N qu10_50_fastGWA_step1
#$ -o /data/scratch/bty207/gcta/step1/qu10_50_logs/
#$ -t 1-69
#$ -tc 13

# GCTA1.94.1 (fastGWA-GLMM) step 1 (GRM) for qu10_50 diseases. Need to partition the GRM into parts for memory reasons (1.35TB required in total: split into 50 parts, requiring max. 30GB each).
# 2 cores requested to speed up GRM generation (would take ~250h with one core).
# Task concurrency set to 25 due to limits of (kindly expanded) scratch storage (have 19.89T free; the full dense GRM takes up ~700G)

mapfile -t diseases < ~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_disease_fields.txt;;
disease=${diseases[${SGE_TASK_ID}]};

module load plink;

# Extract variant set and keep relevant participants, remove chrX. Write the PLINK file set to $TMPDIR
plink2 \
--bfile ~/data/external/genomics/ukb_genotype_data/ukb22418_autosomes-cX_b0_v2 \
--keep ~/data/internal/cohort_eids/genomic_pan-ukbb-eur_eids.tsv \
--extract ~/data/internal/genomics/ld_pruning/${disease}/${disease}_ld_pruning.prune.in \
--not-chr X \
--output-chr MT \
--make-bed \
--threads $NSLOTS \
--out $TMPDIR/autosomes_${disease}-for-grm;

mkdir "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}"; 

# make the full dense GRM, writing the consitutent parts to scratch storage
for i in {1..50};
do ./gcta-1.94.1 \
--bfile $TMPDIR/autosomes_${disease}-for-grm \
--make-grm-part 50 $i \
--threads $NSLOTS \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full";
done;

# remove the PLINK file set from $TMPDIR
rm $TMPDIR/autosomes_${disease}-for-grm*;

# merge parts of full dense GRM, writing the concatenated set to scratch storage ($TMPDIR sometimes runs out of space)
cat "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.part_50_*.grm.N.bin" > "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.grm.N.bin"
cat "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.part_50_*.grm.id" > "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.grm.id"
cat "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.part_50_*.grm.bin" > "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full.grm.bin"

# make sparse GRM
./gcta-1.94.1 \
--grm "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/${disease}_grm_full" \
--make-bK-sparse 0.05 \
--threads $NSLOTS \
--out "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}_grm_sparse"

# delete the concatenated set from $TMPDIR
# rm $TMPDIR/${disease}_grm_full.grm*

mkdir "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/empty_dir_for_deletions/";

# delete the constituent parts of the full dense GRM from scratch storage
time rsync -a --delete "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/empty_dir_for_deletions/" "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}/";
rmdir "~/ch2_genomics/2.2 trial_gwas_runs/2.2.3 fastgwa-glmm/step1/${disease}"
