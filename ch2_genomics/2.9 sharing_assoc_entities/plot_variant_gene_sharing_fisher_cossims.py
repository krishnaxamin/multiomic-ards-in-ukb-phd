"""
Plot variant and gene sharing and cosine similarities on the same set of axes.
Correlations are only calculated for variant-sharing disease pairs as there is no statistic from gene-level MAGMA that
makes sense to correlate.
"""
from pyprind import ProgBar
from typing import List, Tuple, Dict, Any
from utils.significance_labelling import significance_labelling
from utils.get_number_of_variants_tested import get_vars_tested_per_disease_pair, read_in_and_process_information

import pandas as pd
import itertools
import collections
import scipy.stats
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import os
import re
import sklearn.metrics

mpl.use('TkAgg')

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

logger = logging.getLogger('Plotting variant and gene sharing')
logger.setLevel(logging.INFO)
logger.propagate = False

logger.handlers.clear()  # clear existing handlers, including streaming to console

stream_handler = logging.StreamHandler(stream=sys.stdout)
stream_handler.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s - %(name)s. %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)


########################################################################################################################
#
# ENTITY OVERLAP FUNCTIONS
#
########################################################################################################################


def get_entity_overlap(
        data_df: pd.DataFrame, disease_col: str, disease_pairs: List[Tuple[str]], entity_col: str
) -> Dict[str, List[str]]:
    """
    For each disease pair in disease_pairs, identify the set of entities shared between the two.
    :param data_df:
    :param disease_col:
    :param disease_pairs:
    :param entity_col:
    :return:
    """
    shared_entity_dict = collections.defaultdict(list)
    bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Identifying shared entities')
    for disease_pair in disease_pairs:
        disease1 = disease_pair[0]
        disease2 = disease_pair[1]

        shared_entities = list(set(data_df[data_df[disease_col] == disease1][entity_col].unique()) &
                               set(data_df[data_df[disease_col] == disease2][entity_col].unique()))
        shared_entity_dict[disease_pair] = shared_entities
        bar.update()

    return shared_entity_dict


def count_entity_overlap(overlapping_entities: Dict[str, List[str]]) -> pd.DataFrame:
    """
    Produce a DataFrame counting the number of shared entities per disease pair.
    :param overlapping_entities: dictionary from get_entity_overlap, listing shared entities per disease pair.
    :return:
    """
    entity_overlap_df = pd.DataFrame()
    for k, v in overlapping_entities.items():
        entity_overlap_df = pd.concat([entity_overlap_df, pd.DataFrame([{'disease1': k[0],
                                                                         'disease2': k[1],
                                                                         'n_shared': len(v)}])])

    return entity_overlap_df


########################################################################################################################
#
# OVERLAP TESTING
#
########################################################################################################################

def fisher_testing(
        disease_pair: Tuple[str], sharing_counts: pd.DataFrame, entity: str, data_dict,
        n_assoc_entities_per_disease: pd.Series):
    """
    For a specified disease pair, perform Fisher-based overlap testing.
    :param disease_pair:
    :param sharing_counts: DataFrame from count_entity_overlap - the number of shared entities per disease pair
    :param entity: variant or gene
    :param data_dict: if entity==variant: dictionary containing (1) disease info
    (2) dictionary of pd.Series of unique genotyped variant IDs per disease
    (3) pd.Series of unique imputed variant IDs.
    If entity==gene: dictionary of (1) MAGMA results (2) number of total genes testable by MAGMA
    :param n_assoc_entities_per_disease: Series of number of variants/genes associated with each disease
    :return:
    """
    disease1 = disease_pair[0]
    disease2 = disease_pair[1]

    num_shared = int(sharing_counts.loc[(sharing_counts['disease1'] == disease1) &
                                        (sharing_counts['disease2'] == disease2), 'n_shared'].values[0])
    if num_shared == 0:
        return {}

    if entity == 'variant':
        if data_dict is None:
            raise ValueError('data_dict for variant Mainali is not specified.')
        num_tested_both = int(len(get_vars_tested_per_disease_pair(disease1, disease2, data_dict)))

    elif entity == 'gene':
        if data_dict is None:
            raise ValueError('data_dict for gene Mainali is not specified.')
        num_tested_both = int(data_dict['n_total_genes'])

    else:
        raise ValueError(f"Unknown entity: {entity}.")

    # num_blocks_disease1 = len(ld_block_hits[ld_block_hits['disease'] == disease1]['ld_block'].unique())
    # num_blocks_disease2 = len(ld_block_hits[ld_block_hits['disease'] == disease2]['ld_block'].unique())
    #
    # num_shared_blocks = shared_blocks_data_row['count']

    num_assoc_disease1 = int(n_assoc_entities_per_disease[disease1])
    num_assoc_disease2 = int(n_assoc_entities_per_disease[disease2])
    logger.info(f"Fisher inputs: {disease_pair}, {int(num_shared)}, {int(num_assoc_disease1)}, "
                f"{int(num_assoc_disease2)}, {int(num_tested_both)}")

    # set up contingency table [[a, b], [c, d]] where columns are in-test-set +/-, rows are annotation +/-
    #                    disease 1
    #                -------------------
    #                |   +   |   -     |
    #  --------------|-------|---------|
    #   disease 2 +  |   a   |   b     |  num_assoc_disease2
    #   disease 2 -  |   c   |   d     |
    #  --------------|-------|---------|
    #                 num_assoc_disease1  n_tested_both

    a = num_shared
    b = num_assoc_disease2 - num_shared
    c = num_assoc_disease1 - num_shared
    d = num_tested_both - num_assoc_disease1 - b

    fisher_result = scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')

    return {'disease1': disease1, 'disease2': disease2, 'n_shared': num_shared,
            'odds_ratio': fisher_result.statistic, 'pval': fisher_result.pvalue}


def multiple_fisher_testing(disease_pairs: List[Tuple[str]], data_dict, overlap_info: pd.DataFrame, entity: str
                            ) -> pd.DataFrame:
    """
    Perform Fisher-based overlap testing for multiple disease pairs.
    :param disease_pairs:
    :param data_dict: if entity==variant: dictionary containing (1) disease info
    (2) dictionary of pd.Series of unique genotyped variant IDs per disease
    (3) pd.Series of unique imputed variant IDs.
    If entity==gene: dictionary of (1) MAGMA results (2) number of total genes testable by MAGMA
    :param overlap_info: DataFrame from count_entity_overlap - the number of shared entities per disease pair
    :param entity: variant or gene
    :return:
    """
    fisher_dicts = []

    n_assoc_entities_per_disease = data_dict['original_data']['disease'].value_counts()
    bar = ProgBar(len(disease_pairs), stream=sys.stdout, title=f"Fisher testings for '{entity}' workflow")
    for disease_pair in disease_pairs:
        single_fisher_testing = fisher_testing(disease_pair=disease_pair, sharing_counts=overlap_info,
                                               entity=entity, data_dict=data_dict,
                                               n_assoc_entities_per_disease=n_assoc_entities_per_disease)
        if single_fisher_testing:
            fisher_dicts.append(single_fisher_testing)
        bar.update()

    fisher_df = pd.DataFrame(fisher_dicts)

    return fisher_df


########################################################################################################################
#
# CORRELATION CALCULATIONS
#
########################################################################################################################


def get_cossim_of_shared_entities(shared_entities_dicts, overlap_counts: Dict[Any, List[str]],
                                  original_data: pd.DataFrame, original_data_disease_col: str,
                                  original_data_entity_col: str, original_data_effect_col: str) -> pd.DataFrame:
    """
    Given a set of disease pairs which share entities (significance-agnostic), calculate the cosine similarity of those
    entities' effect sizes on each disease.
    :param shared_entities_dicts: DataFrame of disease pairs with entities sharing FDR < 0.05, as a list of dictionaries
    :param overlap_counts: Dict of the entities shared between disease pairs
    :param original_data:
    :param original_data_disease_col: column name in original_data specifying the disease
    :param original_data_entity_col: column name in original_data specifying entity IDs
    :param original_data_effect_col: column name in original_data specifying effect size
    :return:
    """

    data_mat = original_data.pivot(columns=original_data_entity_col, index=original_data_disease_col,
                                   values=original_data_effect_col)

    cossims = pd.DataFrame()
    for shared_entities_dict in shared_entities_dicts:
        disease_pair = (shared_entities_dict['disease1'], shared_entities_dict['disease2'])
        overlapped_entities = overlap_counts[disease_pair]

        overlapped_entities_cossim = sklearn.metrics.pairwise.cosine_similarity(
            data_mat.loc[[shared_entities_dict['disease1'], shared_entities_dict['disease2']], overlapped_entities])

        cossims = pd.concat([cossims,
                             pd.DataFrame([{'disease1': disease_pair[0],
                                            'disease2': disease_pair[1],
                                            'cossim': overlapped_entities_cossim[0, 1]}])])

    return cossims


########################################################################################################################
#
# PLOTTING
#
########################################################################################################################


def preprocess_data_for_heatmaps(input_data: pd.DataFrame, disease_code_type: str, heatmap_triangle: str,
                                 disease_info: pd.DataFrame):
    """
    Organise dataset-to-plot to put the two datasets on opposite halves of the heatmap and recode/recast some variables.
    :param input_data:
    :param disease_code_type:
    :param heatmap_triangle:
    :param disease_info:
    :return:
    """
    # heatmap_triangle = upper, lower. upper = y < x, lower = y > x
    # input_data has disease1 < disease2. heatmap will be plotted with disease1 on x-axis ->
    #  input_data defaults to lower triangle

    processed_data = input_data.copy()

    if heatmap_triangle == 'lower':
        processed_data = processed_data.rename(columns={'disease1': 'disease2', 'disease2': 'disease1'})

    processed_data = (processed_data
                      .merge(disease_info[[disease_code_type, 'code_chapter']]
                             .rename(columns={disease_code_type: 'disease1'}))
                      .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                      .merge(disease_info[[disease_code_type, 'code_chapter']]
                             .rename(columns={disease_code_type: 'disease2'}))
                      .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
    processed_data['disease1'] = pd.Categorical(processed_data['disease1'],
                                                categories=sorted(processed_data['disease1'].unique()))
    processed_data['disease2'] = pd.Categorical(processed_data['disease2'],
                                                categories=sorted(processed_data['disease2'].unique()))
    processed_data['fdr_sig'] = processed_data['fdr_sig'].astype(str)
    processed_data.loc[processed_data['fdr_sig'] == '1', 'fdr_sig'] = 'Yes'
    processed_data.loc[processed_data['fdr_sig'] == '0', 'fdr_sig'] = 'No'

    return processed_data


def plot_combined_heatmap(variant_data_dict, genes_data_dict, plot_path: str, disease_info: pd.DataFrame):
    """
    Plot the Fisher overlap results, with hue=fdr_sig, size=log10(n_shared)
    :param variant_data_dict: Dict containing variant overlap results and other metadata
    :param genes_data_dict: Dict containing gene overlap results and other metadata
    :param plot_path: export path
    :param disease_info:
    :return:
    """
    # showing all instances of overlap, colour=fdr_sig, size=log10(n_shared)

    variant_data_to_plot = preprocess_data_for_heatmaps(variant_data_dict['combined_results'],
                                                        variant_data_dict['disease_code_type'],
                                                        variant_data_dict['heatmap_triangle'],
                                                        disease_info)
    genes_data_to_plot = preprocess_data_for_heatmaps(genes_data_dict['combined_results'],
                                                      genes_data_dict['disease_code_type'],
                                                      genes_data_dict['heatmap_triangle'],
                                                      disease_info)

    data_to_plot = pd.concat([variant_data_to_plot, genes_data_to_plot])

    # deal with cossims
    data_to_plot['cossim'] = data_to_plot['cossim'].fillna(0)
    data_to_plot['cossim_cats'] = data_to_plot['cossim']
    data_to_plot.loc[data_to_plot['cossim'] > 0, 'cossim_cats'] = '~ 1'
    data_to_plot.loc[data_to_plot['cossim'] < 0, 'cossim_cats'] = '~ -1'
    data_to_plot.loc[data_to_plot['cossim'] == 0, 'cossim_cats'] = 'N/A'

    # map data to symmetric categorical index space
    all_labels = sorted(set(data_to_plot['disease1']) | set(data_to_plot['disease2']))
    # make mapping have 1-tick gaps between chapters, so that there are whitespace gaps/boundaries between chapters
    label_to_num = {}
    current_pos = 0
    for i in range(len(all_labels)):
        if i == 0:
            label_to_num[all_labels[i]] = current_pos
            current_pos += 1
            continue
        current_chapter = int(re.search(r"\((\d+)\)", all_labels[i]).group(1))
        past_chapter = int(re.search(r"\((\d+)\)", all_labels[i - 1]).group(1))
        if past_chapter < current_chapter:
            current_pos += 1  # insert a 1-tick wide gap into the axis at the boundary between chapters
        label_to_num[all_labels[i]] = current_pos
        current_pos += 1
    # label_to_num = {lab: i for i, lab in enumerate(all_labels)}
    data_to_plot['disease1_as_num'] = data_to_plot['disease1'].map(label_to_num)
    data_to_plot['disease2_as_num'] = data_to_plot['disease2'].map(label_to_num)

    # plot
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.scatterplot(data=data_to_plot.assign(log_n_shared=np.log10(data_to_plot['n_shared'])).rename(
        columns={'fdr_sig': 'FDR < 0.05', 'log_n_shared': 'log10(number \nvariants/genes \nshared)',
                 'cossim_cats': 'Cosine \nsimilarity'}
    ),
        x='disease1_as_num', y='disease2_as_num', hue='Cosine \nsimilarity', style='FDR < 0.05', style_order=['Yes', 'No'],
        markers={'Yes': 'o', 'No': '^'}, edgecolor=None, palette=["#B40426", "#3B4CC0", 'grey'], hue_order=['~ 1', '~ -1', 'N/A'],
        size='log10(number \nvariants/genes \nshared)', sizes=(25, 160),
        zorder=3, ax=ax)
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
    plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
    plt.grid(True, zorder=0)

    # set ticks labels to the string versions
    tick_positions = [label_to_num[label] for label in all_labels]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(all_labels, rotation=90, size=9)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(all_labels, size=9)

    # add diagonal line going top-left to bottom-right, ensuring whitespaces are dealt with properly
    coords = np.array(tick_positions)
    ax.plot(coords, coords,
            color='black', linestyle='--', linewidth=0.5)

    # bold legend section titles
    legend = ax.get_legend()
    for text, handle in zip(legend.texts, legend.legend_handles):
        if handle._label in ['FDR < 0.05', 'log10(number \nvariants/genes \nshared)', 'Cosine \nsimilarity']:  # section header
            text.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(f"{plot_path}.png")
    plt.savefig(f"{plot_path}.svg")
    plt.close()

    return


########################################################################################################################
#
# ORGANISER AND EXECUTION
#
########################################################################################################################

def data_workflow(input_data: pd.DataFrame, disease_col: str, entity_col: str, effect_col: str, entity: str, **kwargs):
    """
    Use pre-processed association data (variants or genes) to perform Fisher overlap testing and Spearman correlations.
    :param input_data:
    :param disease_col:
    :param entity_col:
    :param effect_col:
    :param entity:
    :param kwargs:
    :return:
    """
    if entity == 'variant':
        data_dict = read_in_and_process_information(logger=logger)
        data_dict['original_data'] = input_data
    elif entity == 'gene':
        all_genes_tested = kwargs.get('all_genes', None)
        if all_genes_tested is None:
            raise ValueError("all_genes tested for MAGMA is not provided.")
        data_dict = {'original_data': input_data, 'n_total_genes': all_genes_tested[5].nunique()}
    else:
        raise ValueError(f"Unknown entity: {entity}.")

    disease_pairs = list(itertools.combinations(sorted(list(input_data[disease_col].unique())), 2))

    overlap = get_entity_overlap(data_df=input_data, disease_col='disease', disease_pairs=disease_pairs,
                                 entity_col=entity_col)
    overlap_counts = count_entity_overlap(overlap)
    logger.info(f"Number of shared entities per disease pair obtained for '{entity}' workflow.")

    fisher = multiple_fisher_testing(disease_pairs=disease_pairs,
                                     data_dict=data_dict,
                                     overlap_info=overlap_counts, entity=entity)
    fisher_sig_labelled = significance_labelling(fisher)
    # fisher_sig_labelled.to_csv(entity_overlap_output_path, index=False)
    logger.info(f"Mainali testing complete for '{entity}' workflow.")

    if entity == 'variant':
        # fisher_sig = fisher_sig_labelled[fisher_sig_labelled['fdr_sig'] == 1].copy()
        fisher_dicts = fisher_sig_labelled.to_dict(orient='records')
        overlapped_cossims = get_cossim_of_shared_entities(shared_entities_dicts=fisher_dicts,
                                                           overlap_counts=overlap,
                                                           original_data=input_data,
                                                           original_data_disease_col=disease_col,
                                                           original_data_entity_col=entity_col,
                                                           original_data_effect_col=effect_col)

        combined_results = fisher_sig_labelled.merge(overlapped_cossims)
    else:  # gene-level MAGMA doesn't have anything to cossim
        combined_results = fisher_sig_labelled.copy()

    # return {'overlap_results': fisher_sig_labelled, 'combined_results': combined_results}
    return combined_results


def execution():
    """
    Main execution function
    :return:
    """
    disease_info = pd.read_csv(
        '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
        str) + ')'

    """ Variants workflow """
    # import if already run
    if os.path.exists(
            '~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_cossims.csv'):
        variants_results_data = pd.read_csv(
            '~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_cossims.csv')
        logger.info('Preloaded variant results.')
    else:

        variants = pd.read_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_lifestyles.csv',
                               dtype={'CHR': str, 'POS': int})
        variants = variants[(variants['assoc_label'] == 'assoc') & (variants['info_label'] == 'high_info')]
        variants['alleles_sorted'] = variants['UKB_ref'] + '_' + variants['UKB_alt']
        variants['alleles_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                      variants['alleles_sorted'].tolist()]
        variants['uniqueID_sorted'] = variants['CHR'].astype(str) + '_' + variants['POS'].astype(
            str) + '_' + variants['alleles_sorted']
        variants.loc[variants['Allele2_is_UKB_alt'] == 0, 'BETA'] = -1 * variants.loc[
            variants['Allele2_is_UKB_alt'] == 0, 'BETA']
        variants = (variants[['disease', 'uniqueID_sorted', 'BETA']]
                    .drop_duplicates(subset=['disease', 'uniqueID_sorted'])
                    .reset_index(drop=True))
        logger.info('Variants read in.')

        variants_results_data = data_workflow(input_data=variants, disease_col='disease', entity_col='uniqueID_sorted',
                                              effect_col='BETA', entity='variant', logger=logger)
        # export
        variants_results_data.to_csv(
            f"~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_cossims.csv",
            index=False)

    """ Genes workflow """
    if os.path.exists(
            '~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher.csv'):
        genes_results_data = pd.read_csv(
            '~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher.csv')
        logger.info('Preloaded gene results.')
    else:
        # Fisher - MAGMA was done using NCBI37.3.gene.loc, so there is a defined universe of genes
        all_genes = pd.read_csv('/data/external/genomics/magma/NCBI37.3/NCBI37.3.gene.loc',
                                sep='\s+', header=None)
        genes_assoc = pd.read_csv('~/data/internal/genomics/magma/lifestyles/gene_level_results_assoc.csv')
        genes_assoc = genes_assoc[genes_assoc['bonf_sig'] == 1][['GENE', 'ZSTAT', 'disease']]
        genes_results_data = data_workflow(input_data=genes_assoc, disease_col='disease', entity_col='GENE',
                                           effect_col='ZSTAT',
                                           entity='gene', all_genes=all_genes)
        # export
        genes_results_data.to_csv(
            f"~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher.csv",
            index=False)

    """ Plotting """
    variants_results_dict = {'combined_results': variants_results_data, 'disease_code_type': 'icd10_three_letter', 'heatmap_triangle': 'upper'}
    genes_results_dict = {'combined_results': genes_results_data, 'disease_code_type': 'disease_field', 'heatmap_triangle': 'lower'}
    plot_combined_heatmap(variant_data_dict=variants_results_dict, genes_data_dict=genes_results_dict,
                          plot_path=f"~/ch2_genomics/2.9 sharing_assoc_entities/plots/"
                                    f"variant_gene_sharing_fisher_testing_cossims",
                          disease_info=disease_info)


if __name__ == '__main__':
    execution()
