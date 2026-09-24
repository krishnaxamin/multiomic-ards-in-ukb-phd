import itertools
import pandas as pd
import numpy as np
import networkx as nx
import scipy
import statsmodels.stats.multitest as smm

association = 'Firth'

""" Louvain """
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
field_to_code_mapper = dict(zip(disease_info.disease_field, disease_info.icd10_three_letter))
code_chapter_to_field_mapper = dict(zip(disease_info.code_chapter, disease_info.disease_field))
bma_similarities = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_disease_resnik_bma.csv")

bma_for_graph = (bma_similarities[['disease1', 'disease2', 'bma_sem_sim']].merge(
    disease_info[['disease_field', 'code_chapter']]
    .rename(columns={'disease_field': 'disease1'}))
                 .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                 .merge(disease_info[['disease_field', 'code_chapter']]
                        .rename(columns={'disease_field': 'disease2'}))
                 .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
bma_graph = nx.Graph()

# [(i, {'disease': <code_chapter>})]
if association == 'firth':
    all_nodes = [(i, {'disease': x}) for i, x in enumerate(disease_info.code_chapter) if x != 'G30 (6)']
else:
    all_nodes = [(i, {'disease': x}) for i, x in enumerate(disease_info.code_chapter)]
bma_graph.add_nodes_from(all_nodes)

disease_node_mapper = {x[1]['disease']: x[0] for x in all_nodes}
node_disease_mapper = {v: k for k, v in disease_node_mapper.items()}

all_edges = []
for i in range(len(bma_for_graph)):
    disease1 = bma_for_graph.at[i, 'disease1']
    disease2 = bma_for_graph.at[i, 'disease2']
    bma_sem_sim = bma_for_graph.at[i, 'bma_sem_sim']

    disease1_node = disease_node_mapper[disease1]
    disease2_node = disease_node_mapper[disease2]

    all_edges.append((disease1_node, disease2_node, {'bma_sem_sim': bma_sem_sim}))

bma_graph.add_edges_from(all_edges)

louvain_clusters = nx.algorithms.community.louvain_communities(bma_graph, seed=42, weight='bma_sem_sim')

# take clusters and Louvain each one
subclusters = []
for cluster in louvain_clusters:
    cluster_graph = bma_graph.subgraph(cluster)
    louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='bma_sem_sim')
    subclusters = subclusters + louvain_subclusters

# more Louvain - break down the clusters as far as they will go
subclusters_next = []
for cluster in subclusters:
    cluster_graph = bma_graph.subgraph(cluster)
    louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='bma_sem_sim')
    subclusters_next += louvain_subclusters

while len(subclusters_next) != len(subclusters):
    subclusters = subclusters_next
    subclusters_next = []
    for cluster in subclusters:
        cluster_graph = bma_graph.subgraph(cluster)
        louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='bma_sem_sim')
        subclusters_next += louvain_subclusters


""" Process Louvain results for fully reduced Louvain clusters """

# isolate non-trivial clusters
louvain_groups = [[node_disease_mapper[x] for x in group] for group in subclusters_next if len(group) > 1]

# is Resnik BMA higher within clusters than without? Valid disease groups should only be those where this is true
louvain_groups_sem_sims = {}
cluster_combinations = []  # disease combinations which are present within clusters
for group in louvain_groups:
    group_diseases_combinations = list(itertools.combinations(group, 2))
    group_sem_sims = []
    for combo in group_diseases_combinations:
        bma_row = bma_for_graph.loc[(bma_for_graph['disease1'] == sorted(combo)[0]) & (
                bma_for_graph['disease2'] == sorted(combo)[1]), 'bma_sem_sim']
        if len(bma_row) > 0:
            group_sem_sims.append(bma_row.values[0])
    louvain_groups_sem_sims['-'.join(group)] = {'sem_sim': group_sem_sims,
                                                'mean_sem_sim': np.mean(group_sem_sims),
                                                'median_sem_sim': np.median(
                                                    group_sem_sims)}  # all medians are above 0.45
    cluster_combinations = cluster_combinations + group_diseases_combinations

cluster_combinations_idxs_bma = []
for combo in cluster_combinations:
    sem_sim_row = bma_for_graph.loc[(bma_for_graph['disease1'] == sorted(combo)[0]) & (
            bma_for_graph['disease2'] == sorted(combo)[1])]
    if len(sem_sim_row) > 0:
        cluster_combinations_idxs_bma.append(sem_sim_row.index.values[0])

non_cluster_sem_sims = bma_for_graph.drop(cluster_combinations_idxs_bma)['bma_sem_sim'].values

louvain_groups_greater_sem_sim = {}
for group in louvain_groups:
    louvain_groups_greater_sem_sim['-'.join(group)] = scipy.stats.mannwhitneyu(
        louvain_groups_sem_sims['-'.join(group)]['sem_sim'],
        non_cluster_sem_sims, alternative='greater').pvalue

# is multimorbidity higher within clusters than without - might this not suggest that multimorbidity is causal for related enrichments?
if association == 'firth':
    multimorbidity = pd.read_csv('~/data/internal/multimorbidity/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences_alphamle')
elif association == 'cox':
    multimorbidity = pd.read_csv('~/data/internal/multimorbidity/proteomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
multimorbidity = pd.merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
    columns={'disease_field': 'disease1'})).rename(columns={'disease1': 'disease1_field', 'code_chapter': 'disease1'})
multimorbidity = pd.merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
    columns={'disease_field': 'disease2'})).rename(columns={'disease2': 'disease2_field', 'code_chapter': 'disease2'})

multimorbidity['disease_pair'] = multimorbidity.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                      axis=1)
multimorbidity = multimorbidity.drop_duplicates('disease_pair')

multimorbidity_fdr_correction = smm.multipletests(multimorbidity['alpha_hat_pval'], alpha=0.05, method='fdr_bh')
multimorbidity = multimorbidity[multimorbidity_fdr_correction[0]].reset_index(drop=True).drop(
    columns=['disease1_field', 'disease2_field'])

multimorbidity[['disease1', 'disease2']] = multimorbidity['disease_pair'].str.split('-', expand=True)

# multimorbidity = multimorbidity.rename(columns={'disease2': 'disease1', 'disease1': 'disease2'})

# get multimorbidity scores for all disease pairs within each disease group
louvain_groups_multimorbidity = {}
# cluster_combinations = []  # disease combinations which are present within clusters
for group in louvain_groups:
    group_diseases_combinations = list(itertools.combinations(group, 2))
    group_alpha_hats = []
    for combo in group_diseases_combinations:
        multimorbidity_combo = multimorbidity.loc[
            (multimorbidity['disease1'] == sorted(combo)[0]) & (
                    multimorbidity['disease2'] == sorted(combo)[1]), 'alpha_hat']
        if len(multimorbidity_combo) == 0:
            continue
        alpha_hat = multimorbidity_combo.values[0]
        group_alpha_hats.append(alpha_hat)
    louvain_groups_multimorbidity['-'.join(group)] = {'alpha_hats': group_alpha_hats,
                                                      'mean_alpha_hat': np.mean(group_alpha_hats),
                                                      'median_alpha_hat': np.median(group_alpha_hats)}
    # cluster_combinations = cluster_combinations + group_diseases_combinations

# mean multimorbidity of whole dataset
mean_all_alpha_hat = multimorbidity['alpha_hat'].mean()

# get alpha-hats for all pairs that are not in the same cluster (non-cluster pairs)
cluster_combinations_idxs_alphamle = []
for combo in cluster_combinations:
    multimorbidity_combo = multimorbidity.loc[
        (multimorbidity['disease1'] == sorted(combo)[0]) & (
                multimorbidity['disease2'] == sorted(combo)[1]), 'alpha_hat']
    if len(multimorbidity_combo) == 0:
        continue
    cluster_combinations_idxs_alphamle.append(multimorbidity_combo.index.values[0])

non_cluster_alpha_hats = multimorbidity.drop(cluster_combinations_idxs_alphamle).drop_duplicates('disease_pair')[
    'alpha_hat'].values

louvain_groups_greater_alphamle = {}
for group in louvain_groups:
    louvain_groups_greater_alphamle['-'.join(group)] = scipy.stats.mannwhitneyu(
        louvain_groups_multimorbidity['-'.join(group)]['alpha_hats'], non_cluster_alpha_hats,
        alternative='greater').pvalue

""" Export Louvain groups """
louvain_groups_as_fields = ['-'.join([code_chapter_to_field_mapper[x] for x in group]) for group in louvain_groups]
louvain_group_df = pd.DataFrame.from_dict(louvain_groups_greater_sem_sim, orient='index')
louvain_group_df.columns = ['sig_greater_resnik_bma']
multimorbidity_df = pd.DataFrame.from_dict(louvain_groups_greater_alphamle, orient='index')
multimorbidity_df.columns = ['sig_greater_multimorbidity']
louvain_group_df = pd.concat([louvain_group_df, multimorbidity_df], axis=1)
louvain_group_df['disease_group'] = louvain_groups_as_fields
louvain_group_df['disease_group'] = ['-'.join(sorted(group)) for group in louvain_groups]
louvain_group_df[['disease_group', 'sig_greater_resnik_bma', 'sig_greater_multimorbidity']].to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_resnik_bma_disease_groups.csv",
    index=False)
