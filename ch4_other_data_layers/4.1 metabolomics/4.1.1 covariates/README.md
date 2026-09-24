# 4.1.1.1 Calculating covariate effects
- Effects calculated using `technical_effects_on_X.R` and `lifestyle_effects_on_X.R`
- Results for each ARD in `technical_effects/` and `lifestyle_effects/`.
- Results collated and plots produced from results using `plotting_technical_effects.py` and `plotting_lifestyle_effects.py`.
- Collated effects results in `technical_effects_on_X.csv` and `lifestyle_effects_on_X.csv`
- Factors affecting protein data only regressed out using `regress_out_effects.py`. Residuals collated and processed into `~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_post_regressions_robust_scaled_data.csv`.
  - Factor effects after each round of regression in above bullet in `./effects_during_regressions.csv`.

# 4.1.1.2. Calculating covariate correlations
- Correlations calculated using `correlation_between_covariates.py`.
- Correlation results in `pearsons_cramers_coeffs_between_covariates_heatmap.csv`.