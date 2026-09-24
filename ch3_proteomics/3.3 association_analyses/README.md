# Code
- `install_coxmos.sh`: Install the R package `Coxmos` for Cox regressions.
- `qu10_50_proteomics_*`: Scripts to perform lifestyle-adjusted Firth or Cox logistic regression on unbalanced proteomics data.
- `process_association_results.py`: Collate results for indiviudal diseases, output from `qu10_50_proteomics_*`, and process by adding significance labels.
- `plot_n_assocs_per_disease_and_protein.py`: Plot number of FDR-sig proteins per ARD, and number of ARDs per associated protein.

# Data absent
- `firth/` and `cox/`: Results for individual diseases, output from `qu10_50_proteomics_*`.