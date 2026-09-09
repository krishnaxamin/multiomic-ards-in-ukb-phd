# Code
- `download_hdl_ref.sh`: Download HDL reference panel.
- `hdl.R`: Collate GWAS summary statistics into an HDL-compatible format and perform HDL.
- `hdl.sh`: Execute HDL for every disease pair.
- `extract_results_from_hdl_output.py`: Read in HDL outputs and concatenate into one csv file.
- `hdl_graph_clustering.py`: Cluster ARDs using FDR-significant HDL correlations.

# Data absent
- `out_files/` and `results/`: Outputs from `hdl.R`. 