# Code
- `qu10_50_metabolomics_lifestyles_firth.R`: Conduct metabolite-disease association analyses using prevalent diagnoses and Firth regressions. Output to `firth_results/`.
- `process_firth_association_results.py`: Collate disease-level association results and process with significance labelling. Outputs `~/data/internal/metabolomics/pan-ukbb-eur_assoc_metabolomics_firth.csv`.
- `plot_n_assocs_per_disease_and_metabolite.py`: Graph to show the number of FDR < 0.05 metabolites per disease, and a heatmap of their effect sizes.

# Data
- `firth_results/`: Contains per-disease association results.