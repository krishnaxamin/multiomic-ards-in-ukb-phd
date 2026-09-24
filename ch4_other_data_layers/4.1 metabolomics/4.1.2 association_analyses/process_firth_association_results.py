""" Process Firth logistic regression metabolomics results. """
from pandas import read_csv, DataFrame, concat
import os
from utils.significance_labelling import significance_labelling


result_files = os.scandir('~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/firth_results')
results_file_names = [x.name for x in result_files if x.is_file()]

results = DataFrame()
for file in results_file_names:
    results = concat([results, read_csv(f"~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/firth_results/{file}")])

results_processed = significance_labelling(results, fdr_permissive_group_by='disease')
results_processed_nom_sig = results_processed[results_processed['nom_sig'] == 1]
results_processed_nom_sig.to_csv(f"~/data/internal/metabolomics/pan-ukbb-eur_assoc_metabolomics_firth.csv", index=False)
