# Code
- `prep_annotation_files.py`: Convert Reactome annotation file annotating Reactome terms to NCBI gene IDs to a MAGMA-compatible format.
- `magma_step3_gene_set_analysis.sh`: Perform MAGMA gene-set analysis using the gene-level results from `magma_step2_collate_summaries_run_magma.sh` and the MAGMA-compatible Reactome annotations from `prep_annotation_files.py`.
- `process_magma_enrichment_results.py`: Collate results from MAGMA enrichment analysis and process for statistical significance. Conduct parentage analysis for select sets of associated pathways.
- `tables_for_disease_pathway_mappings.py`: Generate tables mapping ARDs to MAGMA-associated Reactome terms.