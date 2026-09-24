# Code
- `qu10_50_proteomics_prior_disease_firth.R` and`qu10_50_proteomics_prior_disease_cox.R`: Conduct association analyses that are additionally adjusted for prior disease diagnoses.
- `_master_sub.sh`: Organise HPC jobs to execute association analyses for individual ARD-protein pairs.
- `_nested_subs.sh`: Called by `_master_sub.sh` to execute association analyses for individual ARD-protein pairs.
- `_collate_nested_subs.R`: Collate the individual ARD-protein pairs into results for each ARD.
- `process_association_results.py`: Collate the per-ARD results and process into `pan-ukbb-eur_assoc_proteomics_XXX_prior_disease.csv`.

# Data
- `XXX_results/`: Results from association analyses for individual ARD-protein pairs and collated per-ARD results.