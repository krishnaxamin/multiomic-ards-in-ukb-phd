# Code
- `download_lava_ref.sh`: Download download UK Biobank LD reference. Used in `step5a_preload_lava_input.R`.
- `step0_generate_loci.R`: Generate LD block loci coordinates.
- `step1_sumstats_to_lava.R`: Collate GWAS summary statistics into LAVA-compatible format.
- `step2_same_disease_ldsc.sh`: Conduct same-disease LDSC as required for LAVA. Previous LDSC calculations (in `~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/`) calculated correlations between different diseases only.
- `step3_generate_sample_overlap_file.R`: Generate sample overlap file as required for LAVA from same-disease LDSC results, output from `step2_same_disease_ldsc.sh`.
- `step4_generate_info_file.R`: Generate input info file as required for LAVA, that contains, for each phenotype: the number of cases; the number of controls; the path to the collated summary statistics for that phenotype.
- `step5a_preload_lava_input.R`: Pre-load the LAVA input object to avoid continual loading for each individual array job (LAVA was run on an HPC, with each locus run as a separate job).
- `step5b_run_lava_for_ld_blocks.R`: Perform LAVA for LD blocks.
- `step5b_run_lava_for_follow_up_genes.R`: Perform LAVA for selected genes.
- `collate_lava_results.py`: Collate LAVA heritability and correlation results.
- `analyse_lava_results.py`: Analyse and plot LD block-level LAVA correlation results. Identify genes within specific LD blocks to perform LAVA on. Analyse and plot gene-level LAVA correlation results.

# Data
- `int_data/`: Intermediary data produced and used by the LAVA pipeline.
- `results/`: Locus-level LAVA results for correlation and heritability.