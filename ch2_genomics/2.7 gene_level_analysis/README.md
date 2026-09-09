# Code
- `magma_step1_annotate.sh`: Annotate SNPs from a European reference panel to genes from the NCBI gene set.
- `magma_step2_collate_summaries.R`: Collate GWAS summary statistics into a MAGMA-compatible format.
- `magma_step2_collate_summaries_run_magma.sh`: Run GWAS summary statistics collation, then MAGMA gene-level analysis using the outputs.
- `process_magma_gene_level_results.py`: Collate results from MAGMA gene-level analysis and process for statistical significance.
- `plot_n_diseases_per_gene.py`: Plot how many diseases genes are associated with via MAGMA, with spotlight on those with n_diseases >= 5. Produce table mapping ARDs to MAGMA-associated genes.

# Data absent
- `gwas_summary_stats_for_magma/`: Directory to contain outputs from `magma_step2_collate_summaries.R`.