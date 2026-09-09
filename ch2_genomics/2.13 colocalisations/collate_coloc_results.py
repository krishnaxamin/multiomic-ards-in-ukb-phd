""" Collate ABF and SuSiE colocalisation results for export from Apocrita - run in JupyterLab"""

import os
import pandas as pd

general_path = f"~/ch2_genomics/2.13 colocalisations/results"

for algo in ['abf', 'susie']:
    all_files = [x.name for x in os.scandir(f"{general_path}/{algo}/results") if x.is_file()]
    results = [x for x in all_files if 'results' in x]
    summaries = [x for x in all_files if 'summary' in x]

    pd.concat([pd.read_csv(f"{general_path}/{algo}/results/{x}") for x in results]).to_csv(f"{general_path}/coloc_results_{algo}.csv", index=False)
    pd.concat([pd.read_csv(f"{general_path}/{algo}/results/{x}") for x in summaries]).to_csv(f"{general_path}/coloc_summaries_{algo}.csv", index=False)
