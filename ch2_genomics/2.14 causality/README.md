# Code
- `lcv.R`: Script to run LCV for a particular disease pair.
- `lcv.sh`: Shell script to execute `lcv.R` for all disease pairs.
- `process_lcv.R`: Collate LCV results for each disease pair and label according to statistical significance and 'causality' designation. Output to `~/data/internal/genomics/causality/lcv_results_*`.
- `lcv_icd10_chapter_relationships.py`: Test whether LCV results' relationship to ICD-10 chapters.
- `centrality_metrics.py`: Construct LCV results as a DAG and calculate centrality metrics for each ARD/node from it.
- `lcv_as_graph.dot`: GraphViz dot file to graphically represent the LCV results as a DAG.

# Data absent
- `results/`: Disease pair-level results from `lcv.R`.