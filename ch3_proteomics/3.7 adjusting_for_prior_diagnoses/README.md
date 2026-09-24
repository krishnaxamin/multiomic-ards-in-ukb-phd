# Code
- `get_prior_disease_pcs.py`: Script to isolate prior-disease diagnoses, generate PCs from the diagnoses, and select the number of PCs based on the cumulative proportion of variance explained.
- `get_prior_disease_pcs_proteomics.sh`: Shell script to execute `get_prior_diseases.py` for use in proteomics Firth (prevalent disease) and Cox (incident disease) regressions.
- `prior_disease_pc_XXX_effects_on_phenotype.*`: Script and shell script to calculate the effect of the selected prior-disease PCs on prevalent and incident phenotype.
- `plotting_prior_disease_effects.py`: Plot the calculated effects of prior-disease PCs on phenotype to determine which PCs to include as covariates in the association analyses.
- `association_analyses/`: Contains scripts to perform prior-disease-adjusted association analyses.
- `pathway_based_similarities.py`: Perform Reactome ORA on proteins significantly associated with ARDs through the analyses in `association_analyses/`. Using those enrichment results, calculate pathway-based semantic similarities between ARDs. Cluster the ARDs using the calculated similarities.

# Data
- `XXX_effect_on_phenotype/`: Data of prior-disease PCs' effect on prevalent and incident diagnoses.