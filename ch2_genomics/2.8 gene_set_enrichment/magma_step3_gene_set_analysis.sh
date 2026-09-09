#!bin/bash
#$ -cwd
#$ -l h_rt=1:0:0
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -N magma_step3_gene_sets
#$ -o ~/gwas/magma/magma_step3_logs/
#$ -t 1-68

# MAGMA step 3: gene set analysis

mapfile -t disease_fields < ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt;
disease_field=${disease_fields[${SGE_TASK_ID}]};

# set which GWAS adjustment to use
if [[ "$2" = "minimal" ]] || [[ "$2" = "lifestyles" ]]; then
	adjustment="$2"
else 
	echo "2nd argument must be 'minimal' or 'lifestyles'"
	exit 1
fi

# navigate to directory where the MAGMA executable is
# cd ~/gwas/magma
./magma \
--gene-results ~/data/internal/genomics/magma/${adjustment}/gene_level_results/${disease_field}.genes.raw \
--set-annot ~/data/internal/magma/reactome_annotations_for_magma.tsv col=1,2 \
--model alpha=0.05 \
--out ~/data/internal/genomics/magma/${adjustment}/enrichment_results/${disease_field}
