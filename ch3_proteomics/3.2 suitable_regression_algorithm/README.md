# Code
- `qu10_50_proteomics_smote_rus.py`: Perform SMOTE-RUS balancing and robust scaling on original protein data.
- `qu10_50_proteomics_smote_rus_logistic_regression.R`: Logistic regression performed on SMOTE-RUS-balanced protein data. Standard logistic regression performed for all cases. Elastic net-penalised logistic regression performed only for cases that have _p_ < 0.05, or if a warning arises during normal logistic regression.
- `qu10_50_proteomics_minimal_firth.R`: Framework script that can be altered to run Firth on unbalanced or SMOTE-RUS-balanced data.

# Data absent
- `balanced_data/`: SMOTE-RUS-balanced data for each disease-protein pair, output from `qu10_50_proteomics_smote_rus.py`.
- `association_results/`: Results from normal logistic regression and elastic net-penalised logistic regression, output from `qu10_50_proteomics_smote_rus_logistic_regression.R`.
