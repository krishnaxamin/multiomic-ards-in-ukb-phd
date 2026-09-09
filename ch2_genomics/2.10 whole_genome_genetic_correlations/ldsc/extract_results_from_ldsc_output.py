""" Extract LDSC genetic correlation results from log output """

import glob
import os
import pandas as pd
import argparse
import collections

rows = []

warning_dict = collections.defaultdict(int)
for path in glob.glob("~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/results/*.log"):   # adjust pattern if needed

    output = pd.read_table(path, header=None, delimiter=None)
    output_as_list = list(output[0])

    disease_pair = os.path.basename(path).split('.')[0]
    disease1 = disease_pair.split('-')[0]
    disease2 = disease_pair.split('-')[1]

    warning = [x for x in output_as_list if 'WARNING' in x]
    if len(warning) > 0:
        warning_message = warning[0]
        warning_explainer = output_as_list[output_as_list.index(warning_message) + 1]
        print(f"{disease_pair} flagged '{warning_message} {warning_explainer}")
        warning_dict[warning_message] += 1
        continue

    genetic_correlation = float([x.split(' ')[2] for x in output_as_list if 'Genetic Correlation:' in x][0])
    genetic_correlation_se = float(
        [x.split('(')[1].split(')')[0] for x in output_as_list if 'Genetic Correlation:' in x][0])
    genetic_correlation_pval = float([x.split('P: ')[1].rstrip() for x in output_as_list if 'P: ' in x][0])
    zscore = float([x.split('Z-score: ')[1].rstrip() for x in output_as_list if 'Z-score: ' in x][0])

    rows.append({
        "disease1": disease1,
        "disease2": disease2,
        "genetic_correlation": float(genetic_correlation),
        "se": float(genetic_correlation_se),
        "z": float(zscore),
        "pval": float(genetic_correlation_pval),
    })

for k, v in warning_dict.items():
    print(f"'{k}' had {v} occurrences")

genetic_correlations = pd.DataFrame(rows)
genetic_correlations.to_csv("~/data/internal/genomics/whole_genome_genetic_correlations/ldsc_genetic_correlations.csv", index=False)
