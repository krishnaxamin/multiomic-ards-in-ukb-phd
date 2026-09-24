""" Function to label a dataframe with nominal, FDR, FDR-permissive and Bonferroni significance levels. """
from pandas import DataFrame, concat
from statsmodels.stats.multitest import fdrcorrection


def significance_labelling(input_data, pvalue_col='pval', alpha=0.05, fdr_permissive_group_by=None, pair_mirroring=False,
                           pair_columns=None):
    input_data = input_data.copy()

    # if the input_data contains both permutations of a pair whose data are the same, i.e. not independent tests,
    #  then remove the duplication to accurately capture the number of independent tests
    if not pair_mirroring:
        processing_data = input_data.copy()
    else:
        input_data['pair_sorted'] = input_data.apply(
            lambda row: '_'.join(sorted([row[pair_columns[0]], row[pair_columns[1]]])), axis=1)
        processing_data = input_data.drop_duplicates(['pair_sorted', pvalue_col], ignore_index=True).copy()

    # add significance labels
    processing_data.loc[processing_data[pvalue_col] < alpha, 'nom_sig'] = 1
    processing_data['fdr_sig'], processing_data['pval_fdr_corrected'] = fdrcorrection(processing_data[pvalue_col], alpha=alpha)
    processing_data.loc[processing_data[pvalue_col] < alpha / len(processing_data), 'bonf_sig'] = 1

    # the more generous FDR verion (correcting within groups only, rather than across the whole dataset)
    if fdr_permissive_group_by is not None:
        processing_data_processed = DataFrame()
        for group in processing_data[fdr_permissive_group_by].unique():
            cluster_results = processing_data[processing_data[fdr_permissive_group_by] == group].copy()
            cluster_results['fdr_sig_permissive'], cluster_results['pval_fdr_corrected_permissive'] = fdrcorrection(
                cluster_results[pvalue_col], alpha=alpha)
            processing_data_processed = concat([processing_data_processed, cluster_results])
    else:
        processing_data_processed = processing_data.copy()
    processing_data_processed = processing_data_processed.replace(True, 1).replace(False, 0).fillna(0)

    processing_data_processed['nom_sig'] = processing_data_processed['nom_sig'].astype(int)
    processing_data_processed['bonf_sig'] = processing_data_processed['bonf_sig'].astype(int)

    # if the input_data contains both permutations of a pair whose data are the same, i.e. not independent tests,
    #  restore both permutations of the pair to the output_data
    if not pair_mirroring:
        output_data = processing_data_processed.copy()
    else:
        output_data_med = concat([processing_data_processed,
                              processing_data_processed.assign(disease1=processing_data_processed['disease2'],
                                                               disease2=processing_data_processed['disease1'])])
        output_data = DataFrame()
        for pair_sorted in output_data_med['pair_sorted'].unique():
            output_data = concat([output_data, output_data_med[output_data_med['pair_sorted'] == pair_sorted]])
        output_data.drop(columns='pair_sorted', inplace=True)
        # output_data = processing_data_processed.drop(columns=pair_columns).merge(input_data, on=[x for x in input_data.columns if x not in pair_columns], how='right').drop(columns='pair_sorted')

    return output_data
