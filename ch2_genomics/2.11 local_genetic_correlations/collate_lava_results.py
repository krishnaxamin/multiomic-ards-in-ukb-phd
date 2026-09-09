""" Collate results from individual LAVA runs. """
import pandas as pd
import os

# LD block results
for result_type in ['heritability', 'correlation']:
    files = [x.name for x in os.scandir(f"~/ch2_genomics/2.11 local_genetic_correlations/results/{result_type}") if x.is_file() and 'block' in x.name]
    print(f"Results for LAVA {result_type} analyses number {len(files)}.")
    results = pd.DataFrame()
    for file in files:
        results = pd.concat([results,
                             pd.read_csv(f"~/ch2_genomics/2.11 local_genetic_correlations/results/{result_type}/{file}")])
    results.to_csv(f"~/data/internal/genomics/local_genetic_correlations/lava_ld_block_{result_type}.csv", index=False)

# results for genes selected for follow-up
for result_type in ['heritability', 'correlation']:
    files = [x.name for x in os.scandir(f"~/ch2_genomics/2.11 local_genetic_correlations/results/{result_type}") if x.is_file() and 'ENSG' in x.name]
    print(f"Results for LAVA {result_type} analyses number {len(files)}.")
    results = pd.DataFrame()
    for file in files:
        results = pd.concat([results,
                             pd.read_csv(f"~/ch2_genomics/2.11 local_genetic_correlations/results/{result_type}/{file}")])
    results.to_csv(f"~/data/internal/genomics/local_genetic_correlations/lava_gene_{result_type}.csv", index=False)
