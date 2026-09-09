# Code
- `prep_sumstats_for_ldsc.py`: Collate GWAS summary statistics into LDSC-compatible format.
- `ldsc_prep_sumstats.sh`: Collate GWAS summary statistics into LDSC-compatible format using `prep_sumstats_for_ldsc.py` and munge.
- `ldsc.sh`: Execute LDSC for every disease pair.
- `extract_results_from_ldsc_output.py`: Read in LDSC outputs and concatenate into one csv file.
- `ldsc_graph_clustering.py`: Cluster ARDs using FDR-significant LDSC correlations.

# Data absent
- `gwas_summary_stats_for_ldsc`: Contains subdirs `collated/` and `munged/` respectively containing GWAS summary stats collated ino LDSC-compatible format and munged data to double-check format compatibility.
- `results/`: Outputs from `ldsc.sh`. 
