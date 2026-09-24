from pyprind import ProgBar
from scipy.stats import mannwhitneyu

from utils.significance_labelling import significance_labelling
# from utils.mainali_alpha import mainali_alpha
from utils.enrichment_analyses import overrepresentation, slim_annotations, propagation_up_ontology_hierarchy, \
    prep_annotation_files
from utils.semantic_similarity import semantic_similarity, resnik_bma_between_two_annotation_sets

import pandas as pd
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
import numpy as np
import re
import networkx as nx

mpl.use('TkAgg')

""" Define task """
protein_association = 'firth'

""" Read invariant data """
# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
disease_info_slim = disease_info[['disease_field', 'icd10_three_letter', 'icd10_chapter']].copy()
code_chapter_to_field_mapper = dict(zip(disease_info['code_chapter'], disease_info['disease_field']))

""" Read in ontology data """
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
hierarchy = pd.read_csv(
    '~/data/external/reactome/reactome_id_ancestor_relations.csv').set_axis(
    ['id', 'ancestors'], axis=1)

""" Read in data to use to assess commonalities """
genetic_data = pd.read_csv('~/data/internal/genomics/magma/lifestyles/reactome_enrichment_results_full.csv')
protein_data = pd.read_csv(f"~/data/internal/proteomics/enrichment_analysis/{protein_association}/reactome_ora.csv")
# impute infinite OR to be the maximum finite OR
protein_data.loc[protein_data['odds_ratio'] == np.inf, 'odds_ratio'] = \
    protein_data[protein_data['odds_ratio'] != np.inf]['odds_ratio'].max()

""" Semantic-based similarity """

# calculate IC for all annotations
annotations_ic = -np.log(annotations.annotation_id.value_counts() / len(annotations.id.unique()))

# get absolute value of log10(OR) - used for ranking the top most involved pathways
genetic_data['abs_beta'] = abs(genetic_data['BETA'])
protein_data['abs_log_or'] = abs(np.log10(protein_data['odds_ratio']))

# for each disease group, for each combination of annotations across the diseases, calculate a similarity metric
n_diseases_per_group = 2
# permutation with replacement
disease_pairs = list(
    itertools.product(sorted(list(set(genetic_data.disease))), sorted(list(set(protein_data.disease)))))
genetic_disease_annotations = {disease: d.sort_values('abs_beta', ascending=False, ignore_index=True).loc[:19,
                                        ['VARIABLE']]['VARIABLE'].to_list() for
                               disease, d in genetic_data.groupby('disease')}
protein_disease_annotations = {disease: d.sort_values('abs_log_or', ascending=False, ignore_index=True).loc[:19,
                                        ['annotation']]['annotation'].to_list() for
                               disease, d in protein_data.groupby('disease')}

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

    bma_results, mica_results = resnik_bma_between_two_annotation_sets(annot1=genetic_disease_annotations[disease1],
                                                                       annot2=protein_disease_annotations[disease2],
                                                                       annotation_info_content=annotations_ic,
                                                                       annotation_hierarchy=hierarchy)

    bma_similarities = pd.concat([bma_similarities, bma_results.assign(disease1=disease1, disease2=disease2)])

    # mica_contributions = identify_driver_micas(sem_sim_array=sem_sim_storage_array, mica_array=mica_storage_array)
    bma_mica_contributions = pd.concat([bma_mica_contributions, mica_results.assign(disease1=disease1,
                                                                                    disease2=disease2)])

    bar.update()

""" Label the MICAs with annotation names """
reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                header=None, names=['id', 'name', 'species'])
reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
bma_mica_contributions_export = bma_mica_contributions.merge(
    reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

""" Export semantic similarities """
bma_similarities.to_csv(
    f"~/data/internal/integration/gwas_{protein_association}_enrichment_resnik_bma.csv",
    index=False)
bma_mica_contributions_export.to_csv(
    f"~/data/internal/integration/gwas_{protein_association}_enrichment_resnik_micas.csv",
    index=False)

""" Heatmap """
pairwise_similarities = bma_similarities[['disease1', 'disease2', 'bma_sem_sim']].copy()
pairwise_similarities_to_plot = (pairwise_similarities
                                 .merge(disease_info[['disease_field', 'code_chapter']]
                                        .rename(columns={'disease_field': 'disease1'}))
                                 .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                                 .merge(disease_info[['disease_field', 'code_chapter']]
                                        .rename(columns={'disease_field': 'disease2'}))
                                 .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
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
plt.xlabel('Disease ICD-10 code (ICD-10 chapter) - Genetic Enrichments')
plt.ylabel('Disease ICD-10 code (ICD-10 chapter) - Proteomic (' + protein_association[0].upper() + protein_association[
                                                                                                   1:] + ') Enrichments')
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
plt.savefig(
    f"~/ch5_integration/5.1 pathway_based_similarities/plots/heatmap_disease_pairs_gwas_{protein_association}_enrichment_resnik_bma.png")
plt.savefig(
    f"~/ch5_integration/5.1 pathway_based_similarities/plots/heatmap_disease_pairs_gwas_{protein_association}_enrichment_resnik_bma.svg")
plt.close()

""" Analysis """
# expect better alignment from same-disease pairs than diff-disease pairs
mannwhitneyu(
    bma_similarities[bma_similarities['disease1'] == bma_similarities['disease2']][
        'bma_sem_sim'],
    bma_similarities[bma_similarities['disease1'] != bma_similarities['disease2']][
        'bma_sem_sim'])

# same-disease semsim from Firth-GWAS vs Cox-GWAS
firth_bma = pd.read_csv(f"~/data/internal/integration/gwas_firth_enrichment_resnik_bma.csv")
cox_bma = pd.read_csv(f"~/data/internal/integration/gwas_cox_enrichment_resnik_bma.csv")
mannwhitneyu(cox_bma[cox_bma['disease1'] == cox_bma['disease2']]['bma_sem_sim'],
             firth_bma[firth_bma['disease1'] == firth_bma['disease2']]['bma_sem_sim'], alternative='greater')
