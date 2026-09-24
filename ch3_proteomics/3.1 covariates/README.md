# 3.1.1. Calculating covariate effects
- Effects calculated using `technical_effects_on_X.R` and `lifestyle_effects_on_X.R`
- Results for each ARD in `technical_effects/` and `lifestyle_effects/`.
- Collated effects results in `technical_effects_on_X.csv` and `lifestyle_effects_on_X.csv`
- Plots produced from results using `plotting_technical_effects.py` and `plotting_lifestyle_effects.py`.
- Factors affecting protein data only regressed out using `regress_out_effects.py` and `regress_out_effects_for_cox.py`.
  - Residuals from above bullet in `regressing_out/post_regression_protein_data/` and `regressing_out_for_cox/post_regression_protein_data/`.
  - Factor effects after each round of regression in above bullet in `regressing_out/effects/` and `regressing_out_for_cox/effects/`.
  - Collated factor effects in `regressing_out/effects_during_regressions.csv` and `regressing_out_for_cox/effects_during_regressions_for_cox.csv`. Data absent to reduce storage requirements.

# 3.1.2. Calculating covariate correlations
- Correlations calculated using `correlation_between_covariates.py`.
- Correlation results in `pearsons_cramers_coeffs_between_covariates_heatmap.csv`.