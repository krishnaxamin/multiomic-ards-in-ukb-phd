"""
Perform graph-based Louvain clustering FDR-corrected significant LDSC values.
Iterate Louvain clustering on clusters until no new clusters are formed.
Discard any trivial networks (single nodes) that result.
==> disease groups
Calculate significance of correlations and multimorbidity within clusters being greater than non-cluster relationships.
Export disease groups.
"""
from pandas import read_csv, concat, DataFrame, merge
from itertools import combinations
import statsmodels.stats.multitest as smm
from scipy.stats import mannwhitneyu

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx

mpl.use('TkAgg')

ldsc = read_csv(f"~/data/internal/genomics/whole_genome_genetic_correlations/ldsc_genetic_correlations.csv")
ldsc.loc[ldsc['genetic_correlation'] > 1, 'genetic_correlation'] = 1
ldsc.loc[ldsc['genetic_correlation'] < -1, 'genetic_correlation'] = -1
disease_info = read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

# dictionary for switching between disease field and disease ICD-10 code and chapter, and vice versa
field_to_code_chapter_mapper = dict(zip(disease_info['disease_field'], disease_info['code_chapter']))
code_chapter_to_field_mapper = dict(zip(disease_info['code_chapter'], disease_info['disease_field']))

# FDR correction and filtering
fdr_correction = smm.multipletests(ldsc['pval'], alpha=0.05, method='fdr_bh')
ldsc_sig = ldsc[fdr_correction[0]]

""" Set up graph """

ldsc_for_graph = merge(ldsc_sig,
                       disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease1',
                                                                                       'code_chapter': 'disease1_code_chapter'}),
                       on='disease1')
ldsc_for_graph = merge(ldsc_for_graph,
                       disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease2',
                                                                                       'code_chapter': 'disease2_code_chapter'}),
                       on='disease2')

all_graph = nx.Graph()

# [(i, {'disease': <code_chapter>})]
all_nodes = [(i, {'disease': x}) for i, x in enumerate(disease_info.code_chapter)]
all_graph.add_nodes_from(all_nodes)

node_disease_mapper = dict(zip(list(range(68)), list(disease_info.code_chapter)))
disease_node_mapper = dict(zip(list(disease_info.code_chapter), list(range(68))))

all_edges = []
for i in range(len(ldsc_for_graph)):
    disease1 = ldsc_for_graph.at[i, 'disease1_code_chapter']
    disease2 = ldsc_for_graph.at[i, 'disease2_code_chapter']
    cor = ldsc_for_graph.at[i, 'genetic_correlation']

    disease1_node = disease_node_mapper[disease1]
    disease2_node = disease_node_mapper[disease2]

    all_edges.append((disease1_node, disease2_node, {'cor': cor, 'abs_cor': abs(cor)}))

all_graph.add_edges_from(all_edges)

""" Louvain """  # incorporates edge weights
louvain_clusters = nx.algorithms.community.louvain_communities(all_graph, seed=42, weight='cor')

# take clusters and Louvain each one
subclusters = []
for cluster in louvain_clusters:
    cluster_graph = all_graph.subgraph(cluster)
    louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='cor')
    subclusters = subclusters + louvain_subclusters

# more Louvain - break down the clusters as far as they will go
subclusters2 = []
for cluster in subclusters:
    cluster_graph = all_graph.subgraph(cluster)
    louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='cor')
    subclusters2 = subclusters2 + louvain_subclusters
subclusters3 = []
for cluster in subclusters2:
    cluster_graph = all_graph.subgraph(cluster)
    louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='cor')
    subclusters3 = subclusters3 + louvain_subclusters

# subclusters3 == subclusters2 so that's the max number of Louvain-based clusters

disease_groups = [[node_disease_mapper[x] for x in group] for group in [list(x) for x in subclusters3] if len(group) > 1]

disease_groups_correlations = {}
cluster_combinations = []  # disease combinations which are present within clusters
for group in disease_groups:
    group_diseases_combinations = list(combinations(group, 2))
    group_corrs = []
    for combo in group_diseases_combinations:
        corr_row = ldsc_for_graph.loc[(ldsc_for_graph['disease1_code_chapter'] == sorted(combo)[0]) & (
                    ldsc_for_graph['disease2_code_chapter'] == sorted(combo)[1]), 'genetic_correlation']
        if len(corr_row) > 0:
            group_corrs.append(corr_row.values[0])
    disease_groups_correlations['-'.join(group)] = {'corrs': group_corrs,
                                                    'mean_corr': np.mean(group_corrs),
                                                    'median_corr': np.median(group_corrs)}  # all medians are above 0.45
    cluster_combinations = cluster_combinations + group_diseases_combinations

cluster_combinations_idxs_ldsc = []
for combo in cluster_combinations:
    corr_row = ldsc_for_graph.loc[(ldsc_for_graph['disease1_code_chapter'] == sorted(combo)[0]) & (
                ldsc_for_graph['disease2_code_chapter'] == sorted(combo)[1])]
    if len(corr_row) > 0:
        cluster_combinations_idxs_ldsc.append(corr_row.index.values[0])

non_cluster_corrs = ldsc_for_graph.drop(cluster_combinations_idxs_ldsc)['genetic_correlation'].values

disease_groups_greater_corr = {}
for group in disease_groups:
    disease_groups_greater_corr['-'.join(group)] = mannwhitneyu(disease_groups_correlations['-'.join(group)]['corrs'],
                                                                non_cluster_corrs, alternative='greater').pvalue

""" Are multimorbidities within the final cohorts significantly stronger than non-cluster multimorbidities? """
multimorbidity = read_csv('~/data/internal/multimorbidity/genomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
multimorbidity_fdr_correction = smm.multipletests(multimorbidity['alpha_hat_pval'], alpha=0.05, method='fdr_bh')
multimorbidity = multimorbidity[multimorbidity_fdr_correction[0]].reset_index(drop=True)

multimorbidity = merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
    columns={'disease_field': 'disease1'})).rename(columns={'disease1': 'disease1_field', 'code_chapter': 'disease1'})
multimorbidity = merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
    columns={'disease_field': 'disease2'})).rename(columns={'disease2': 'disease2_field', 'code_chapter': 'disease2'})

multimorbidity['disease_pair'] = multimorbidity.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                      axis=1)

# get multimorbidity scores for all disease pairs within each disease group
disease_groups_multimorbidity = {}
# cluster_combinations = []  # disease combinations which are present within clusters
for group in disease_groups:
    group_diseases_combinations = list(combinations(group, 2))
    group_alpha_hats = []
    for combo in group_diseases_combinations:
        alpha_hat = multimorbidity.loc[
            (multimorbidity['disease1'] == combo[0]) & (multimorbidity['disease2'] == combo[1]), 'alpha_hat'].values[0]
        group_alpha_hats.append(alpha_hat)
    disease_groups_multimorbidity['-'.join(group)] = {'alpha_hats': group_alpha_hats,
                                                      'mean_alpha_hat': np.mean(group_alpha_hats),
                                                      'median_alpha_hat': np.median(group_alpha_hats)}
    # cluster_combinations = cluster_combinations + group_diseases_combinations

# mean multimorbidity of whole dataset = 1.252
mean_all_alpha_hat = multimorbidity['alpha_hat'].mean()

# get alpha-hats for all pairs that are not in the same cluster (non-cluster pairs)
cluster_combinations_idxs_alphamle = []
for combo in cluster_combinations:
    cluster_combinations_idxs_alphamle.append(multimorbidity.loc[(multimorbidity['disease1'] == combo[0]) & (
                multimorbidity['disease2'] == combo[1])].index.values[0])
    cluster_combinations_idxs_alphamle.append(multimorbidity.loc[(multimorbidity['disease2'] == combo[0]) & (
                multimorbidity['disease1'] == combo[1])].index.values[0])

non_cluster_alpha_hats = multimorbidity.drop(cluster_combinations_idxs_alphamle).drop_duplicates('disease_pair')[
    'alpha_hat'].values

disease_groups_greater_alphamle = {}
for group in disease_groups:
    disease_groups_greater_alphamle['-'.join(group)] = mannwhitneyu(
        disease_groups_multimorbidity['-'.join(group)]['alpha_hats'], non_cluster_alpha_hats,
        alternative='greater').pvalue

""" Export disease groups with correlation and multimorbidity info """
disease_groups_as_str = ['-'.join(sorted(group)) for group in disease_groups]
corr_df = DataFrame.from_dict(disease_groups_greater_corr, orient='index')
corr_df.columns = ['sig_greater_correlation']
multimorbidity_df = DataFrame.from_dict(disease_groups_greater_alphamle, orient='index')
multimorbidity_df.columns = ['sig_greater_multimorbidity']
disease_group_df = concat([corr_df, multimorbidity_df], axis=1)
disease_group_df['disease_group'] = disease_groups_as_str
disease_group_df[['disease_group', 'sig_greater_correlation', 'sig_greater_multimorbidity']].to_csv(
    "~/data/internal/genomics/whole_genome_genetic_correlations/ldsc_disease_groups.csv", index=False)
