from pyprind import ProgBar

from utils.enrichment_analyses import prep_annotation_files
from utils.semantic_similarity import resnik_bma_between_two_annotation_sets

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy
import numpy as np
import networkx as nx
import pandas as pd
import itertools
import statsmodels.stats.multitest as smm
import math
import re

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

""" Semantic-based similarity """

# calculate IC for all annotations
annotations_ic = -np.log(annotations.annotation_id.value_counts() / len(annotations.id.unique()))

# calculate weights for all annotation-disease pairs
data['weights'] = -np.log10(data['P']) * np.log10(data['BETA'])

# get absolute value of log10(OR) - used for ranking the top most involved pathways
data['abs_beta'] = abs(data['BETA'])

# for each disease group, for each combination of annotations across the diseases, calculate a similarity metric
disease_pairs = list(itertools.combinations(sorted(list(data['disease'].unique())), 2))
disease_annotation_sets = {disease: d.sort_values('abs_beta', ascending=False, ignore_index=True).loc[:19,
                                    ['VARIABLE']]['VARIABLE'].to_list() for disease, d in data.groupby('disease')}

# bma similarities
# for each disease pair A,B, set up df with columns=A terms, index=B terms, populate df with sem-sims
# sem-sims are symmetric: sim(a,b) = sim(b,a)
# column-wise and index-wide max() to get all max sem-sims for each of A and B terms
# sum maxes, divide each sum by number of A and B terms, sum results and divide by 2
# BMA only defined for pairwise similarity

bma_similarities = pd.DataFrame()
bma_mica_contributions = pd.DataFrame()
bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Semantic similarity calculation')
for disease_pair in disease_pairs:
    disease1 = disease_pair[0]
    disease2 = disease_pair[1]

    bma_results, mica_results = resnik_bma_between_two_annotation_sets(annot1=disease_annotation_sets[disease1],
                                                                       annot2=disease_annotation_sets[disease2],
                                                                       annotation_info_content=annotations_ic,
                                                                       annotation_hierarchy=hierarchy)

    bma_similarities = pd.concat([bma_similarities, bma_results.assign(disease1=disease1, disease2=disease2)])

    # mica_contributions = identify_driver_micas(sem_sim_array=sem_sim_storage_array, mica_array=mica_storage_array)
    bma_mica_contributions = pd.concat([bma_mica_contributions, mica_results.assign(disease1=disease1,
                                                                                    disease2=disease2)])

    bar.update()

""" Export semantic similarities """
bma_similarities.to_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_bma_reactome.csv", index=False)
bma_mica_contributions.to_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_micas_reactome.csv",
                              index=False)


""" Within-chapter vs between-chapter SemSims """


def semsims_within_between_chapters_mwu(semsims):
    """
    Perform a Mann-Whitney U rank sum test for the difference in semantic similarity between same-chapter diseases and
    between different-chapter diseases.
    :param semsims:
    :return:
    """
    bma_within_between = (semsims[['disease1', 'disease2', 'bma_sem_sim']]
    .merge(disease_info[['disease_field', 'icd10_chapter']].rename(
        columns={'disease_field': 'disease1', 'icd10_chapter': 'disease1_chapter'}))
    .merge(disease_info[['disease_field', 'icd10_chapter']].rename(
        columns={'disease_field': 'disease2', 'icd10_chapter': 'disease2_chapter'})))
    bma_within_between['chapter_match'] = bma_within_between['disease1_chapter'] == bma_within_between[
        'disease2_chapter']

    return scipy.stats.mannwhitneyu(bma_within_between[bma_within_between['chapter_match']].bma_sem_sim,
                                    bma_within_between[~bma_within_between['chapter_match']].bma_sem_sim,
                                    alternative='greater')


bma_within_between_u_rw = semsims_within_between_chapters_mwu(bma_similarities).statistic

# generate null distribution of Mann-Whitney U for within-vs-between chapter SemSims
bma_within_between_u_null = []
# bma_similarities_for_perms = bma_similarities.copy()
bma_sem_sims = bma_similarities.bma_sem_sim.to_list()
rng = np.random.default_rng(42)
n_iterations = 10000
bar = ProgBar(n_iterations, stream=sys.stdout, title='Generating null for SemSim within-vs-between chapters')
for _ in range(n_iterations):
    shuffled_sem_sims = rng.permutation(bma_sem_sims)
    bma_within_between_u_null.append(
        semsims_within_between_chapters_mwu(bma_similarities.assign(bma_sem_sim=shuffled_sem_sims)).statistic)
    bar.update()

# when comparing RW to null, don't need to normalise to AUC (done by dividing U stat by n1xn2) because n1,n2 always identical through every permutation
# higher U = within-SemSims > between-SemSims
print(
    f"RW >= null empirical p: {1 - sum([x <= bma_within_between_u_rw for x in bma_within_between_u_null]) / n_iterations}")
print(
    f"RW <= null empirical p: {1 - sum([x >= bma_within_between_u_rw for x in bma_within_between_u_null]) / n_iterations}")

""" Heatmap """
pairwise_similarities = bma_similarities[['disease1', 'disease2', 'bma_sem_sim']].copy()
pairwise_similarities_to_plot = (pairwise_similarities
                                 .merge(disease_info[['disease_field', 'code_chapter']]
                                        .rename(columns={'disease_field': 'disease1'}))
                                 .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                                 .merge(disease_info[['disease_field', 'code_chapter']]
                                        .rename(columns={'disease_field': 'disease2'}))
                                 .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
pairwise_similarities_to_plot = pd.concat([pairwise_similarities_to_plot,
                                           pairwise_similarities_to_plot.assign(
                                               disease1=pairwise_similarities_to_plot.disease2,
                                               disease2=pairwise_similarities_to_plot.disease1)])
pairwise_similarities_to_plot['disease1'] = pd.Categorical(pairwise_similarities_to_plot['disease1'],
                                                           categories=sorted(pairwise_similarities_to_plot[
                                                                                 'disease1'].unique()))
pairwise_similarities_to_plot['disease2'] = pd.Categorical(pairwise_similarities_to_plot['disease2'],
                                                           categories=sorted(pairwise_similarities_to_plot[
                                                                                 'disease2'].unique()))

# map data to symmetric categorical index space
all_labels = sorted(set(pairwise_similarities_to_plot['disease1']) | set(pairwise_similarities_to_plot['disease2']))
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
pairwise_similarities_to_plot['disease1_as_num'] = pairwise_similarities_to_plot['disease1'].map(label_to_num)
pairwise_similarities_to_plot['disease2_as_num'] = pairwise_similarities_to_plot['disease2'].map(label_to_num)

fig, ax = plt.subplots(figsize=(12, 10))
sns.scatterplot(pairwise_similarities_to_plot, x='disease1_as_num', y='disease2_as_num', hue='bma_sem_sim',
                palette='viridis', zorder=3, ax=ax)
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', title='Resnik BMA', fontsize=8)
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.grid(True, zorder=0)

if ax.yaxis_inverted():
    ax.invert_yaxis()

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

# build colorbar
norm = mpl.colors.Normalize(vmin=pairwise_similarities_to_plot.bma_sem_sim.min(),
                            vmax=pairwise_similarities_to_plot.bma_sem_sim.max())
sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
sm.set_array([])  # required boilerplate

# remove legend, attach colorbar
ax.get_legend().remove()
plt.colorbar(sm, ax=ax, label='Semantic similarity')

plt.tight_layout()
plt.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/heatmap_disease_pairs_by_reactome_resnik_bma.svg")
plt.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/heatmap_disease_pairs_by_reactome_resnik_bma.png")
plt.close()
