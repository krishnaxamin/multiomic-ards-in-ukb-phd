from pyprind import ProgBar
from scipy.stats import mannwhitneyu
from concepts import Context

from utils.enrichment_analyses import prep_annotation_files
from utils.semantic_similarity import resnik_bma_between_two_annotation_sets

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
import numpy as np
import networkx as nx
import pandas as pd
import itertools
import scipy
import pylluvial as pa
import statsmodels.stats.multitest as smm
import math
import re
import itertools
import upsetplot

mpl.use('TkAgg')

""" Read invariant data """
# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
disease_info_slim = disease_info[['disease_field', 'icd10_three_letter', 'icd10_chapter']].copy()
code_chapter_to_field_mapper = dict(zip(disease_info['code_chapter'], disease_info['disease_field']))
field_to_code_chapter_mapper = dict(zip(disease_info['disease_field'], disease_info['code_chapter']))

""" Read in ontology data """
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/NCBI2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
hierarchy = pd.read_csv(
    '~/data/external/reactome/reactome_id_ancestor_relations.csv').set_axis(
    ['id', 'ancestors'], axis=1)

""" Read in data to use to assess commonalities """
data = pd.read_csv('~/data/internal/genomics/magma/lifestyles/reactome_enrichment_results_full.csv')

# calculate IC for all annotations
annotations_ic = -np.log(annotations.annotation_id.value_counts() / len(annotations.id.unique()))

# calculate weights for all annotation-disease pairs
data['weights'] = -np.log10(data['P']) * np.log10(data['BETA'])

# get absolute value of log10(OR) - used for ranking the top most involved pathways
data['abs_beta'] = abs(data['BETA'])

""" Semantic-based similarity - top-10 """

# for each disease group, for each combination of annotations across the diseases, calculate a similarity metric
disease_pairs = list(itertools.combinations(sorted(list(data['disease'].unique())), 2))
disease_annotation_sets = {disease: d.sort_values('abs_beta', ascending=False, ignore_index=True).loc[:9,
                                    ['VARIABLE']]['VARIABLE'].to_list() for disease, d in data.groupby('disease')}

# bma similarities
# for each disease pair A,B, set up df with columns=A terms, index=B terms, populate df with sem-sims
# sem-sims are symmetric: sim(a,b) = sim(b,a)
# column-wise and index-wide max() to get all max sem-sims for each of A and B terms
# sum maxes, divide each sum by number of A and B terms, sum results and divide by 2
# BMA only defined for pairwise similarity

bma_similarities_top10 = pd.DataFrame()
bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Semantic similarity calculation')
for disease_pair in disease_pairs:
    disease1 = disease_pair[0]
    disease2 = disease_pair[1]

    bma_results, mica_results = resnik_bma_between_two_annotation_sets(annot1=disease_annotation_sets[disease1],
                                                                       annot2=disease_annotation_sets[disease2],
                                                                       annotation_info_content=annotations_ic,
                                                                       annotation_hierarchy=hierarchy)

    bma_similarities_top10 = pd.concat(
        [bma_similarities_top10, bma_results.assign(disease1=disease1, disease2=disease2)])

    bar.update()
bma_similarities_top10.to_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_bma_reactome_top10.csv",
                              index=False)

""" Semantic-based similarity - top-15 """

disease_annotation_sets = {disease: d.sort_values('abs_beta', ascending=False, ignore_index=True).loc[:14,
                                    ['VARIABLE']]['VARIABLE'].to_list() for disease, d in data.groupby('disease')}

# bma similarities
# for each disease pair A,B, set up df with columns=A terms, index=B terms, populate df with sem-sims
# sem-sims are symmetric: sim(a,b) = sim(b,a)
# column-wise and index-wide max() to get all max sem-sims for each of A and B terms
# sum maxes, divide each sum by number of A and B terms, sum results and divide by 2
# BMA only defined for pairwise similarity

bma_similarities_top15 = pd.DataFrame()
bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Semantic similarity calculation')
for disease_pair in disease_pairs:
    disease1 = disease_pair[0]
    disease2 = disease_pair[1]

    bma_results, mica_results = resnik_bma_between_two_annotation_sets(annot1=disease_annotation_sets[disease1],
                                                                       annot2=disease_annotation_sets[disease2],
                                                                       annotation_info_content=annotations_ic,
                                                                       annotation_hierarchy=hierarchy)

    bma_similarities_top15 = pd.concat(
        [bma_similarities_top15, bma_results.assign(disease1=disease1, disease2=disease2)])

    bar.update()
bma_similarities_top15.to_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_bma_reactome_top15.csv",
                              index=False)

""" Semantic similarity - top 20 """
bma_similarities_top20 = pd.read_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_bma_reactome.csv")

""" SemSim comparisons """
bma_similarities_comparison = (
    bma_similarities_top20[['disease1', 'disease2', 'bma_sem_sim']].rename(columns={'bma_sem_sim': 'top20_semsim'})
    .merge(
        bma_similarities_top10[['disease1', 'disease2', 'bma_sem_sim']].rename(columns={'bma_sem_sim': 'top10_semsim'}))
    .merge(bma_similarities_top15[['disease1', 'disease2', 'bma_sem_sim']].rename(
        columns={'bma_sem_sim': 'top15_semsim'})))

# top15 vs top20 correlation
scipy.stats.pearsonr(bma_similarities_comparison['top15_semsim'], bma_similarities_comparison['top20_semsim'])
# top10 vs top20 correlation
scipy.stats.pearsonr(bma_similarities_comparison['top10_semsim'], bma_similarities_comparison['top20_semsim'])

""" Compare Louvain results """


def louvain(bma_similarities):
    # bma_semsim_strong = bma_similarities[bma_similarities['bma_sem_sim'] > bma_similarities['bma_sem_sim'].median()]
    bma_for_graph = (bma_similarities[['disease1', 'disease2', 'bma_sem_sim']].merge(
        disease_info[['disease_field', 'code_chapter']]
        .rename(columns={'disease_field': 'disease1'}))
                     .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                     .merge(disease_info[['disease_field', 'code_chapter']]
                            .rename(columns={'disease_field': 'disease2'}))
                     .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
    bma_graph = nx.Graph()

    # [(i, {'disease': <code_chapter>})]
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
            louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42,
                                                                              weight='bma_sem_sim')
            subclusters_next += louvain_subclusters

    # isolate non-trivial clusters
    recursive_clusters = [[node_disease_mapper[x] for x in group] for group in
                          [list(x) for x in subclusters_next if len(x) > 1]]

    return recursive_clusters


top10_louvain = louvain(bma_similarities_top10)
top15_louvain = louvain(bma_similarities_top15)
top20_louvain = louvain(bma_similarities_top20)

""" Alluvials to plot changes in groupings """
# top10/15/20 = 'timepoint', points on the x axis (doesn't have to be called 'timepoint'). arg: x =
# alluvium = ..., the specific entities
# stratum = ..., the categories the entities fall in at each 'timepoint'
top10_pylluvial = pd.concat(
    [pd.DataFrame({'cat': [i + 1] * len(lst), 'disease': lst}) for i, lst in enumerate(top10_louvain)]).assign(
    timepoint='top10')
top15_pylluvial = pd.concat(
    [pd.DataFrame({'cat': [i + 1] * len(lst), 'disease': lst}) for i, lst in enumerate(top15_louvain)]).assign(
    timepoint='top15')
top20_pylluvial = pd.concat(
    [pd.DataFrame({'cat': [i + 1] * len(lst), 'disease': lst}) for i, lst in enumerate(top20_louvain)]).assign(
    timepoint='top20')

pylluvial_data = pd.concat([top10_pylluvial, top15_pylluvial, top20_pylluvial])
pylluvial_data['timepoint'] = pd.Categorical(pylluvial_data['timepoint'], ['top10', 'top20', 'top15'])

fig, ax = pa.alluvial(
    data=pylluvial_data,
    x='timepoint',
    stratum='cat',
    alluvium='disease',
    show_labels=True
)

# removing edge lines (by setting width = 0) removes the white lines that represent 0 flow between 2 categories
for patch in ax.patches:
    patch.set_linewidth(0)

# expand vertically to ensure capture of all flows
ymin, ymax = ax.get_ylim()
ax.set_ylim(ymin, ymax * 1.05)

fig.set_figwidth(12)
fig.set_figheight(8)
fig.tight_layout()

fig.savefig(f"~/data/internal/genomics/pathway_based_similarities/plots/semsim_diff_topN_louvain_groups.png")
fig.savefig(f"~/data/internal/genomics/pathway_based_similarities/plots/semsim_diff_topN_louvain_groups.svg")
plt.close()

""" FCA to compare Louvain groups """
# Build incidence matrix: rows = groups, columns = diseases
# diseases = sorted({d for g in groups for d in g})
# df = pd.DataFrame(False, index=range(len(groups)), columns=diseases)
# for i, g in enumerate(groups):
#     df.loc[i, list(g)] = True
fca_prep_df = pylluvial_data.copy()
fca_prep_df['group'] = fca_prep_df['timepoint'].astype(str) + '-' + fca_prep_df['cat'].astype(str)
fca_mat_df = fca_prep_df[['group', 'disease']].assign(values=True).reset_index(drop=True).pivot(
    columns='disease', index='group', values='values').fillna(False)

# Run FCA
ctx = Context(objects=fca_mat_df.index.astype(str).tolist(),
              properties=fca_mat_df.columns.astype(str).tolist(),
              bools=fca_mat_df.astype(bool).values.tolist())
lattice = ctx.lattice

# Extract interpretable overlap patterns
overlap_patterns_list = []
bar = ProgBar(len(lattice), stream=sys.stdout,
              title='FCA for disease groups shared between Louvain clusters of different top-N')
for concept in lattice:
    diseases_in_common = concept.intent  # exact disease set
    groups_sharing = concept.extent  # groups that share them

    # Filter trivial concepts
    if len(diseases_in_common) >= 2 and len(groups_sharing) >= 2:
        overlap_patterns_list.append({
            'n_diseases': len(diseases_in_common),
            'n_clusters': len(groups_sharing),
            'diseases': '|'.join(sorted(list(diseases_in_common))),
            'clusters': '|'.join(sorted(list(groups_sharing)))
        })
        # print(f"Diseases {set(diseases_in_common)}"
        #       f"appear in groups {set(groups_sharing)}")

    bar.update()
overlap_patterns_df = pd.DataFrame(overlap_patterns_list)
overlap_patterns_df = overlap_patterns_df.sort_values(by='n_clusters', ascending=False)
overlap_patterns_df.to_csv(f"~/data/internal/genomics/pathway_based_similarities/semsim_diff_topN_louvain_groups_fca_results.csv",
                           index=False)

# number of diseases grouped with at least one other disease over different top-N clusters
len(set(itertools.chain.from_iterable([x.split('|') for x in overlap_patterns_df[overlap_patterns_df.n_clusters == 3].diseases.to_list()])))
overlap_patterns_df.clusters.str.contains('top20').sum()
# number of diseases grouped with at least one other disease in a top20 cluster and another top-N cluster
len(set(itertools.chain.from_iterable([x.split('|') for x in overlap_patterns_df[(overlap_patterns_df.n_clusters == 2) & (overlap_patterns_df.clusters.str.contains('top20'))].diseases.to_list()])))

""" Vary top-N for MICAs"""


def count_high_flying_micas(selected_micas_df: pd.DataFrame) -> pd.DataFrame:
    high_flying_micas_slim = selected_micas_df[['mica', 'disease1', 'disease2']].drop_duplicates()
    high_flying_micas_diseases = [
        {'mica': mica,
         'diseases': '|'.join([field_to_code_chapter_mapper[x]
                               for x in
                               sorted(list(set(mica_diseases_df['disease1']) | set(mica_diseases_df['disease2'])))])}
        for mica, mica_diseases_df in high_flying_micas_slim.groupby('mica')]
    high_flying_micas_diseases_df = pd.DataFrame(high_flying_micas_diseases)
    high_flying_micas_diseases_df['n_diseases'] = high_flying_micas_diseases_df['diseases'].apply(
        lambda x: len(x.split('|')))
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        pd.DataFrame(high_flying_micas_slim['mica'].value_counts().rename('pair_count')).reset_index())

    # get reactome pathway names
    reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                    header=None, names=['id', 'name', 'species'])
    reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

    return high_flying_micas_diseases_df.sort_values(by='pair_count', ascending=False).reset_index(drop=True)


def get_complete_connectors(micas_counted):
    micas_counted_working = micas_counted.copy()
    micas_counted_working['pairwise_completion'] = [
        micas_counted_working.at[i, 'pair_count'] / math.comb(
            micas_counted_working.at[i, 'n_diseases'],
            2) for i in
        range(len(micas_counted_working))]
    micas_counted_working_complete_connectors = micas_counted_working[
        micas_counted_working['pairwise_completion'] == 1].drop(columns='pairwise_completion').copy()
    micas_counted_working_complete_connectors = micas_counted_working_complete_connectors[
        micas_counted_working_complete_connectors['n_diseases'] > 2].copy()

    return micas_counted_working_complete_connectors


bma_mica_contributions_top20 = pd.read_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_micas_reactome.csv")

top5_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:5, :] for _, disease_pair_df in bma_mica_contributions_top20.groupby(['disease1', 'disease2'])])
prioritised_micas_top5 = top5_micas_per_disease_pair[
    top5_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
prioritised_micas_top5_counted = count_high_flying_micas(prioritised_micas_top5)
complete_connectors_top5 = get_complete_connectors(prioritised_micas_top5_counted)

top3_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:3, :] for _, disease_pair_df in bma_mica_contributions_top20.groupby(['disease1', 'disease2'])])
prioritised_micas_top3 = top3_micas_per_disease_pair[
    top3_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
prioritised_micas_top3_counted = count_high_flying_micas(prioritised_micas_top3)
complete_connectors_top3 = get_complete_connectors(prioritised_micas_top3_counted)

top10_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:10, :] for _, disease_pair_df in bma_mica_contributions_top20.groupby(['disease1', 'disease2'])])
prioritised_micas_top10 = top10_micas_per_disease_pair[
    top10_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
prioritised_micas_top10_counted = count_high_flying_micas(prioritised_micas_top10)
complete_connectors_top10 = get_complete_connectors(prioritised_micas_top10_counted)

# upset for numbers of complete connectors
complete_connectors_upset = upsetplot.from_contents({'top3': complete_connectors_top3.mica.to_list(),
                                                     'top5': complete_connectors_top5.mica.to_list(),
                                                     'top10': complete_connectors_top10.mica.to_list()})
fig = plt.figure(figsize=(8, 6))
complete_connectors_upsetplot = upsetplot.UpSet(complete_connectors_upset, element_size=None, show_counts=True,
                                                sort_categories_by='input')
complete_connectors_upsetplot.style_subsets(min_degree=3, facecolor='#1f77b4')  # colour by 1st tab10 colour
complete_connectors_upsetplot.plot(fig=fig)
fig.savefig(f"~/data/internal/genomics/pathway_based_similarities/plots/complete_connectors_diff_topN_micas.png")
fig.savefig(f"~/data/internal/genomics/pathway_based_similarities/plots/complete_connectors_diff_topN_micas.svg")
plt.close()

# n_diseases of complete connectors shared between all 3 topNs
complete_connectors_compare_n_diseases = (complete_connectors_top3[['mica_name', 'n_diseases']].rename(columns={'n_diseases': 'top3'})
                                          .merge(complete_connectors_top5[['mica_name', 'n_diseases']].rename(columns={'n_diseases': 'top5'}))
                                          .merge(complete_connectors_top10[['mica_name', 'n_diseases']].rename(columns={'n_diseases': 'top10'})))