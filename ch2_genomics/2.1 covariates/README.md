# 2.1.1. Calculating genetic principal components
- `generate_genetic_principal_components.sh`

# 2.1.2. Calculating covariate effects
- Effects calculated using `technical_effects_on_phenotype.R` and `lifestyle_effects_on_phenotype.R`
- Results for each ARD in `technical_effects/` and `lifestyle_effects/`.
- Collated effects results in `technical_effects_on_phenotype.csv` and `lifestyle_effects_on_phenotype.csv`
- Plots produced from results using `plotting_technical_effects.py` and `plotting_lifestyle_effects.py`.

# 2.1.3. Calculating covariate correlations
- Correlations calculated using `correlation_between_covariates.py`.
- Correlation results in `pearsons_cramers_coeffs_between_covariates_heatmap.csv`.