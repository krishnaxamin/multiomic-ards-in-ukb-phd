""" Read in HDL outputs and concatenate into one csv file. """
from pandas import DataFrame, read_table, concat
from pyprind import ProgBar

import sys
import os


# for all output files
output_files = [x.name for x in os.scandir('~/ch2_genomics/2.10 whole_genome_genetic_correlations/hdl/results') if x.is_file()]
result_df_list = []
bar = ProgBar(len(output_files), stream=sys.stdout, title='Processing HDL outputs')
for output_file in output_files:
    disease_pair = output_file.split('.')[0]
    disease1 = disease_pair.split('-')[0]
    disease2 = disease_pair.split('-')[1]

    hdl_output = read_table('~/ch2_genomics/2.10 whole_genome_genetic_correlations/hdl/results/' + output_file, header=None, delimiter=None)
    hdl_output_as_list = list(hdl_output[0])
    genetic_correlation = float([x.split(' ')[3] for x in hdl_output_as_list if 'Genetic Correlation' in x][0])
    genetic_correlation_sd = float(
        [x.split('(')[1].split(')')[0] for x in hdl_output_as_list if 'Genetic Correlation' in x][0])
    genetic_correlation_pval = float([x.split('P: ')[1].rstrip() for x in hdl_output_as_list if 'P: ' in x][0])
    result_df = DataFrame([{'disease1': disease1,
                            'disease2': disease2,
                            'genetic_correlation': genetic_correlation,
                            'sd': genetic_correlation_sd,
                            'pval': genetic_correlation_pval}])
    result_df_list.append(result_df)
    bar.update()

results_df = concat(result_df_list).sort_values(by=['disease1', 'disease2'])
results_df.to_csv('~/data/internal/genomics/whole_genome_genetic_correlations/hdl_genetic_correlations.csv', index=False)
