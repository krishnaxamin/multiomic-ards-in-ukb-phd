""" Process lifestyle-and-prior-disease-adjusted regression proteomics results. """
import os

from pandas import read_csv, DataFrame, concat
from utils.significance_labelling import significance_labelling

association = 'firth'  # or 'cox'

result_files = os.scandir(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/{association}_results")
results_file_names = [x.name for x in result_files if x.is_file()]

results = DataFrame()
for file in results_file_names:
    results = concat([results, read_csv(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/{association}_results/{file}")])

results_processed = significance_labelling(results, fdr_permissive_group_by='disease')
results_processed_nom_sig = results_processed[results_processed['nom_sig'] == 1]
results_processed_nom_sig.to_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{association}_prior_disease.csv", index=False)
