""" Quantify similarities between Firth and Cox ORA top-20-by-effect Reactome pathways. Explore certain relationships. """

# looking for high similarity between Firth for one disease and Cox for another
# expect high similarity between Firth and Cox for the same disease?

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
import numpy as np
import pandas as pd
import itertools
import re
import scipy
import math

# from collections import defaultdict
from pyprind import ProgBar
# from scipy.spatial.distance import squareform
# from scipy.cluster.hierarchy import linkage
# from sklearn.preprocessing import MinMaxScaler

from utils.significance_labelling import significance_labelling
# from utils.mainali_alpha import mainali_alpha
from utils.enrichment_analyses import prep_annotation_files, reactome_parentage_analysis
from utils.semantic_similarity import resnik_bma_between_two_annotation_sets

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
field_to_code_mapper = dict(zip(disease_info.disease_field, disease_info.icd10_three_letter))

""" Read in ontology data """
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
hierarchy = pd.read_csv(
    '~/data/external/reactome/reactome_id_ancestor_relations.csv').set_axis(
    ['id', 'ancestors'], axis=1)

""" Read in data to use to assess commonalities """
firth_data = pd.read_csv('~/data/internal/proteomics/enrichment_analysis/firth/reactome_ora.csv')
# impute infinite OR to be the maximum finite OR
firth_data.loc[firth_data['odds_ratio'] == np.inf, 'odds_ratio'] = firth_data[firth_data['odds_ratio'] != np.inf][
    'odds_ratio'].max()

cox_data = pd.read_csv('~/data/internal/proteomics/enrichment_analysis/cox/reactome_ora.csv')
# impute infinite OR to be the maximum finite OR
cox_data.loc[cox_data['odds_ratio'] == np.inf, 'odds_ratio'] = cox_data[cox_data['odds_ratio'] != np.inf][
    'odds_ratio'].max()

""" Read in mapping data """
olink_uniprot_mappings = pd.read_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv')

""" Define protein and annotation background """
# --- protein background ---
# use the mapping file, as this contains all the UKB proteins before QC
# remove the 3 proteins calculated as high missingness during QC
high_missingness_proteins = pd.read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/proteomics/high_missingness_protein_fields.txt')
all_qc_proteins = olink_uniprot_mappings[
    ~olink_uniprot_mappings['protein_field'].isin(high_missingness_proteins['protein_field'])].copy()
protein_background = all_qc_proteins['uniprot'].to_list()

# --- annotation background ---
annotations_for_protein_background = annotations[annotations['id'].isin(protein_background)].copy()
annotations_for_protein_background = annotations_for_protein_background[['id', 'annotation_id']].drop_duplicates()

# slim
# annotations_for_protein_background is only needed to calculate IC
# no need to slim because
#  (1) IC is dependent on genes, not terms, in the data. Slimming affects terms only, note genes
#  (2) slimming would remove high-level MICAs necessary for SemSim calcs
# annotations_for_protein_background = slim_annotations(annotations=annotations_for_protein_background,
#                                                       annotation_col_id='annotation_id', entity_col_id='id')

# calculate IC for all annotations FOR ONLY THOSE PROTEINS ANALYSED
annotations_ic = -np.log(annotations_for_protein_background.annotation_id.value_counts() / len(
    annotations_for_protein_background.id.unique()))

# no need to check that slimming has worked, since the background annotation set hasn't been slimmed
# if wanting to check that slimming is OK, filter to the pathways present in the enrichment results
assert annotations_ic[firth_data.annotation].max() < -np.log(10 / 19026)
assert annotations_ic[cox_data.annotation].max() < -np.log(10 / 19026)

""" Semantic-based similarity """

firth_data['abs_log_or'] = abs(np.log10(firth_data['odds_ratio']))
cox_data['abs_log_or'] = abs(np.log10(cox_data['odds_ratio']))

# for each disease group, for each combination of annotations across the diseases, calculate a similarity metric
# genetic_protein_diseases = set(firth_data.disease) | set(cox_data.disease)
# permutation with replacement
disease_pairs = list(itertools.product(sorted(list(set(firth_data.disease))), sorted(list(set(cox_data.disease)))))
firth_disease_annotation_sets = {disease: d.sort_values('abs_log_or', ascending=False, ignore_index=True).loc[:19,
                                          ['annotation']]['annotation'].to_list() for
                                 disease, d in firth_data.groupby('disease')}
cox_disease_annotation_sets = {disease: d.sort_values('abs_log_or', ascending=False, ignore_index=True).loc[:19,
                                        ['annotation']]['annotation'].to_list() for
                               disease, d in cox_data.groupby('disease')}
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

    bma_results, mica_results = resnik_bma_between_two_annotation_sets(annot1=firth_disease_annotation_sets[disease1],
                                                                       annot2=cox_disease_annotation_sets[disease2],
                                                                       annotation_info_content=annotations_ic,
                                                                       annotation_hierarchy=hierarchy)

    bma_similarities = pd.concat([bma_similarities, bma_results.assign(disease1=disease1, disease2=disease2)])

    # mica_contributions = identify_driver_micas(sem_sim_array=sem_sim_storage_array, mica_array=mica_storage_array)
    bma_mica_contributions = pd.concat([bma_mica_contributions, mica_results.assign(disease1=disease1,
                                                                                    disease2=disease2)])

    bar.update()

""" Label the MICAs with annotation names """
# bma_mica_contributions = pd.read_csv(f"ukbiobank/proteomics/proteomic_similarities/{association}/disease_resnik_micas_{ontology}_{adjustment}.csv")
reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                header=None, names=['id', 'name', 'species'])
reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
bma_mica_contributions_export = bma_mica_contributions.merge(
    reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

""" Export semantic similarities """
bma_similarities.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_disease_resnik_bma.csv",
    index=False)
bma_mica_contributions.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_disease_resnik_micas.csv",
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
                palette='viridis', zorder=3, ax=ax, s=50)
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', title='Resnik BMA', fontsize=8)
plt.xlabel('Prevalent (Firth) disease ICD-10 code (ICD-10 chapter)', size=11)
plt.ylabel('Incident (Cox) disease ICD-10 code (ICD-10 chapter)', size=11)
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
    f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_heatmap_disease_pairs_by_resnik_bma.svg")
plt.savefig(
    f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_heatmap_disease_pairs_by_resnik_bma.png")
plt.close()

""" Within-chapter vs between-chapter SemSims """


def semsims_within_between_chapters_mwu(semsims, alternative='greater'):
    """
    Perform a Mann-Whitney U rank sum test for the difference in semantic similarity between same-chapter diseases and
    between different-chapter diseases.
    :param alternative:
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
                                    alternative=alternative)


semsims_within_between_chapters_mwu(bma_similarities, alternative='two-sided')

""" Same-disease vs different-disease SemSims"""


def semsims_within_between_diseases_mwu(semsims, alternative='greater'):
    """
    Perform a Mann-Whitney U rank sum test for the difference in Firth-Cox semantic similarity between the same diseases
    and between different diseases.
    :param alternative:
    :param semsims:
    :return:
    """
    bma_disease_match = semsims.copy()
    bma_disease_match['disease_match'] = bma_disease_match.disease1 == bma_disease_match.disease2

    return scipy.stats.mannwhitneyu(bma_disease_match[bma_disease_match['disease_match']].bma_sem_sim,
                                    bma_disease_match[~bma_disease_match['disease_match']].bma_sem_sim,
                                    alternative=alternative)


semsims_within_between_diseases_mwu(bma_similarities, alternative='greater')

""" Identify top-contributing MICAs """


def count_high_flying_micas(selected_micas_df: pd.DataFrame) -> pd.DataFrame:
    high_flying_micas_slim = selected_micas_df[['mica', 'disease1', 'disease2']].drop_duplicates()
    high_flying_micas_diseases = [
        {'mica': mica,
         'diseases_firth_to_cox': '|'.join(
             ['-'.join(
                 [field_to_code_chapter_mapper[firth_cox_tup[0]],
                  field_to_code_chapter_mapper[firth_cox_tup[1]]])
                 for firth_cox_tup in list(zip(mica_diseases_df['disease1'], mica_diseases_df['disease2']))]),
         'n_diseases': len(set(mica_diseases_df['disease1']) | set(mica_diseases_df['disease2'])),
         # 'diseases_firth_to_cox': '|'.join([field_to_code_chapter_mapper[x]
         #                       for x in
         #                       sorted(list(set(mica_diseases_df['disease1']) | set(mica_diseases_df['disease2'])))])
         }
        for mica, mica_diseases_df in high_flying_micas_slim.groupby('mica')]
    high_flying_micas_diseases_df = pd.DataFrame(high_flying_micas_diseases)
    # high_flying_micas_diseases_df['n_diseases'] = high_flying_micas_diseases_df['diseases'].apply(
    #     lambda x: len(x.split('|')))
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        pd.DataFrame(high_flying_micas_slim['mica'].value_counts().rename('pair_count')).reset_index())

    # get reactome pathway names
    reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                    header=None, names=['id', 'name', 'species'])
    reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

    return high_flying_micas_diseases_df.sort_values(by='pair_count', ascending=False).reset_index(drop=True)


def plot_high_flying_micas(counted_micas_df: pd.DataFrame, n_top_pathways: int, plot_path: str):
    """
    Given a set of high-flying driver MICAs, plot their links with diseases and disease pairs via a
    (1) heatmap (MICA-disease)
    (2) bar chart counting n(diseases)
    (3) bar chart counting n(disease pairs)
    The MICAs are ordered by n(disease pairs) because we are interested in MICAs that drive similarity between pairs.
    :param counted_micas_df:
    :param n_top_pathways:
    :param plot_name:
    :return:
    """
    # heatmap. mica_name on y-axis (cut-off after n characters), ARDs x-axis.
    # Separate bar charts of n_diseases and pair_count, as with the combined manhattan-ish plot
    counted_micas_df = counted_micas_df.reset_index(drop=True)

    # convert counted_micas_df to heatmap format
    counted_micas_heatmap = pd.DataFrame(index=counted_micas_df['mica_name'], columns=disease_info['code_chapter'],
                                         dtype=int)

    for _, row in counted_micas_df.iterrows():
        counted_micas_heatmap.loc[row['mica_name'], row['diseases_firth_to_cox'].split('|')] = 1
    counted_micas_heatmap.fillna(0, inplace=True)

    top_n_pathways_counts = counted_micas_df.iloc[:n_top_pathways, :].copy()

    # Figure layout
    fig = plt.figure(figsize=(16, 12))
    fig.set_constrained_layout(True)
    gs = fig.add_gridspec(1, 3, width_ratios=[4, 1, 1], wspace=0.05)

    # --- Scatter plot ---
    ax = fig.add_subplot(gs[0])

    sns.heatmap(
        data=counted_micas_heatmap.iloc[:n_top_pathways, :], linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
        cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
        linecolor='black'
    )

    ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
    ax.set_ylabel('Reactome term')
    ax.set_xticklabels(ax.get_xticklabels(), size=8)
    ax.set_yticklabels(ax.get_yticklabels(), size=8)

    # --- Bar chart (counts per disease) ---
    # log_vars = np.log1p(n_vars_per_disease_full.values)
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    ax2.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['n_diseases'], color='grey',
             alpha=0.7)
    ax2.set_xlabel('n(diseases)')
    ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

    # Add vertical grid lines at tick marks
    xticks = ax2.get_xticks()
    for xt in xticks:
        ax2.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['n_diseases'].to_list()):
        ax2.text(count, y, str(count), va='center', fontsize=8)

    # --- Bar chart (gene counts per disease) ---
    # log_genes = np.log1p(n_genes_per_disease_full.values)
    ax3 = fig.add_subplot(gs[2], sharey=ax)
    ax3.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['pair_count'], color='grey',
             alpha=0.7)
    ax3.set_xlabel('n(disease pairs)')
    ax3.tick_params(axis='y', left=False, labelleft=False)

    # Vertical grid lines
    for xt in ax3.get_xticks():
        ax3.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['pair_count'].to_list()):
        ax3.text(count, y, str(count), va='center', fontsize=8)

    plt.show()

    plt.savefig(f"{plot_path}.png")
    plt.savefig(f"{plot_path}.svg")

    plt.close()

    return


# segregate bma_mica_contributions into (1) same-disease (2) Firth-Cox (3) Cox-Firth
bma_mica_contributions_seg = {
    'same_disease': bma_mica_contributions[bma_mica_contributions.disease1 == bma_mica_contributions.disease2].copy(),
    'firth_cox': bma_mica_contributions[bma_mica_contributions.disease1 < bma_mica_contributions.disease2].copy(),
    'cox_firth': bma_mica_contributions[bma_mica_contributions.disease1 > bma_mica_contributions.disease2].copy(),
}

# bma_mica_contributions_between_disease = bma_mica_contributions[bma_mica_contributions.disease1 != bma_mica_contributions.disease2].copy()
# for contribution_label, contribution_df in bma_mica_contributions_seg.items():
# top 5 per disease pair by contribution + MICA IC >= 5.587705 (25th percentile of all Reactome ICs)
top5_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:5, :] for _, disease_pair_df in bma_mica_contributions.groupby(['disease1', 'disease2'])])
prioritised_micas = top5_micas_per_disease_pair[
    top5_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
prioritised_micas_counted = count_high_flying_micas(prioritised_micas)

# look at these MICAs but index by how many disease-disease connections they're involved in given the set of diseases they're linked to
#  so if a MICA is linked to 10 diseases, there's a maximum of 45 pairs - how many of these is the MICA representing/linked to
prioritised_micas_counted['pairwise_completion'] = [
    prioritised_micas_counted.at[i, 'pair_count'] / math.comb(prioritised_micas_counted.at[i, 'n_diseases'],
                                                               2)
    for i in range(len(prioritised_micas_counted))]
(prioritised_micas_counted.pairwise_completion == 1).sum()
prioritised_micas_complete_connectors = prioritised_micas_counted[
    prioritised_micas_counted['pairwise_completion'] == 1].drop(columns='pairwise_completion').copy()
prioritised_micas_complete_connectors = prioritised_micas_complete_connectors[
    prioritised_micas_complete_connectors['n_diseases'] > 2].copy()
prioritised_micas_counted.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_semsim_common_driver_pathways.csv",
    index=False)
prioritised_micas_complete_connectors.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_semsim_common_driver_pathways_complete_connectors.csv",
    index=False)

# plot
if not prioritised_micas_complete_connectors.empty:
    # continue
    plot_path = f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_semsim_common_driver_pathways_complete_connectors"
    plot_high_flying_micas(counted_micas_df=prioritised_micas_complete_connectors, n_top_pathways=25,
                           plot_path=plot_path)

    # parentage analysis
    _ = reactome_parentage_analysis(annotation_set=prioritised_micas_complete_connectors.mica.to_list(),
                                    plot_path=f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_semsim_common_driver_pathways_complete_connectors_parentage",
                                    export_plot=True)

# parentage analysis
_ = reactome_parentage_analysis(annotation_set=prioritised_micas_counted.mica.to_list(),
                                plot_path=f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_semsim_common_driver_pathways_parentage",
                                export_plot=True)

""" Correlation between Firth and Cox Resniks """
firth_resniks = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_disease_resnik_bma.csv")
cox_resniks = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/cox_disease_resnik_bma.csv")

firth_cox_resniks = (
    firth_resniks[['disease1', 'disease2', 'bma_sem_sim']].rename(columns={'bma_sem_sim': 'firth_semsim'})
    .merge(cox_resniks[['disease1', 'disease2', 'bma_sem_sim']].rename(columns={'bma_sem_sim': 'cox_semsim'})))
scipy.stats.spearmanr(firth_cox_resniks.firth_semsim, firth_cox_resniks.cox_semsim)  # 0.122, p = 0.00122
scipy.stats.pearsonr(firth_cox_resniks.firth_semsim, firth_cox_resniks.cox_semsim)  # 0.142, p = 1.55e-4

# plot
firth_cox_resniks_to_plot = (firth_cox_resniks
                             .merge(disease_info[['disease_field', 'code_chapter']]
                                    .rename(columns={'disease_field': 'disease1'}))
                             .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                             .merge(disease_info[['disease_field', 'code_chapter']]
                                    .rename(columns={'disease_field': 'disease2'}))
                             .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
plt.figure(figsize=(8, 6))
sns.scatterplot(firth_cox_resniks, x='firth_semsim', y='cox_semsim')
plt.axline((0, 0), slope=1, linestyle='--', linewidth=0.5, color='black')
plt.xlabel('Firth-Firth semantic similarity')
plt.ylabel('Cox-Cox semantic similarity')
plt.tight_layout()
plt.savefig(f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_correlation.svg")
plt.savefig(f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_cox_correlation.png")
plt.close()

""" Firth ARDs similar to many Cox ARDs """
bma_similarities = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_disease_resnik_bma.csv")
ards_by_chapter = {chapter: df.disease_field.to_list() for chapter, df in disease_info.groupby('icd10_chapter')}
one_firth_many_cox_mw_list = []
for disease in bma_similarities.disease1.unique():
    # for each Firth ARD, test whether its SemSims with Cox ARDs are greater than everyone else
    results_dict = {'firth_disease': disease, 'all_ards': scipy.stats.mannwhitneyu(
        bma_similarities[bma_similarities['disease1'] == disease].bma_sem_sim,
        bma_similarities[bma_similarities['disease1'] != disease].bma_sem_sim,
        alternative='greater').pvalue}
    # same test, but resolved for each chapter, i.e. for each Firth ARD, are they particularly SemSim with Cox ARDs of a particular chapter
    for chapter, ards in ards_by_chapter.items():
        results_dict[f"chapter_{chapter}"] = scipy.stats.mannwhitneyu(
            bma_similarities[
                (bma_similarities['disease1'] == disease) & (bma_similarities['disease2'].isin(ards))].bma_sem_sim,
            bma_similarities[
                (bma_similarities['disease1'] != disease) & (bma_similarities['disease2'].isin(ards))].bma_sem_sim,
            alternative='greater').pvalue
    one_firth_many_cox_mw_list.append(results_dict)

one_firth_many_cox_mw = pd.DataFrame(one_firth_many_cox_mw_list)

one_firth_many_cox_mw_long = one_firth_many_cox_mw.melt(id_vars='firth_disease', var_name='cox_ard_set',
                                                        value_name='greater_pval')

# plot
one_firth_many_cox_mw_to_plot = one_firth_many_cox_mw.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'firth_disease'})).drop(
    'firth_disease', axis=1)
one_firth_many_cox_mw_to_plot.rename(columns={'all_ards': 'All ARDs'} |
                                             {col: f"Chapter {col.split('_')[-1]}"
                                              for col in one_firth_many_cox_mw_to_plot.columns if 'chapter_' in col},
                                     inplace=True)
one_firth_many_cox_mw_to_plot_long = one_firth_many_cox_mw_to_plot.melt(id_vars='code_chapter', var_name='cox_ard_set',
                                                                        value_name='greater_pval')
one_firth_many_cox_mw_to_plot_long = significance_labelling(one_firth_many_cox_mw_to_plot_long,
                                                            pvalue_col='greater_pval')
one_firth_many_cox_mw_to_plot = one_firth_many_cox_mw_to_plot_long[one_firth_many_cox_mw_to_plot_long.fdr_sig == 1][
    ['code_chapter', 'cox_ard_set', 'fdr_sig']].pivot(index='code_chapter',
                                                      columns='cox_ard_set',
                                                      values='fdr_sig').fillna(0)
one_firth_many_cox_mw_to_plot.columns = ['All ARDs', 'Chapter 4', 'Chapter 6', 'Chapter 7', 'Chapter 9', 'Chapter 10',
                                         'Chapter 11', 'Chapter 13']

# Figure layout
fig = plt.figure(figsize=(10, 8))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 7], wspace=0.05)

# --- All ARDs heatmap ---
ax = fig.add_subplot(gs[0])

sns.heatmap(
    data=one_firth_many_cox_mw_to_plot.loc[:, ['All ARDs']], linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
    cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
    vmin=0, vmax=1,
    linecolor='black'
)

ax.set_xlabel('')
ax.set_ylabel('Prevalent disease ICD-10 code (ICD-10 chapter)')
# ax.set_xticklabels(ax.get_xticklabels(), size=7)
ax.set_yticklabels(ax.get_yticklabels(), size=9)

# --- ARDs-by-chapter heatmap ---
ax2 = fig.add_subplot(gs[1], sharey=ax)

sns.heatmap(
    data=one_firth_many_cox_mw_to_plot.drop(columns='All ARDs'), linewidths=0.35, ax=ax2, xticklabels=True,
    yticklabels=True,
    cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
    linecolor='black'
)

ax2.set_xlabel('Incident ARD set')
ax2.set_ylabel('')
# ax2.set_xticklabels(ax2.get_xticklabels(), size=7, rotation=90)
ax2.tick_params(axis='y', left=False, labelleft=False)

plt.savefig(f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_ard_semsim_to_many_cox_ards.png")
plt.savefig(f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firth_ard_semsim_to_many_cox_ards.svg")
plt.close()

""" Chapter 13 one-Firth-many-Cox investigation """
bma_similarities = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/firth_cox_disease_resnik_bma.csv")

# I25, K76, N18
prevalent_ard_set = ['p131306', 'p131670', 'p132032']

# hip and knee OA, vs other OA
scipy.stats.mannwhitneyu(
    bma_similarities[(bma_similarities.disease1.isin(prevalent_ard_set)) & (
        bma_similarities.disease2.isin(['p131870', 'p131872']))].bma_sem_sim,
    bma_similarities[(bma_similarities.disease1.isin(prevalent_ard_set)) & (
        bma_similarities.disease2.isin(['p131868', 'p131874']))].bma_sem_sim,
)

# OAs, vs other chapter 13s
scipy.stats.mannwhitneyu(
    bma_similarities[(bma_similarities.disease1.isin(prevalent_ard_set)) & (
        bma_similarities.disease2.isin(['p131868', 'p131870', 'p131872', 'p131874']))].bma_sem_sim,
    bma_similarities[(bma_similarities.disease1.isin(prevalent_ard_set)) & (
        bma_similarities.disease2.isin(['p131900', 'p131964', 'p131972']))].bma_sem_sim,
    alternative='greater'
)

# albation test - does the signficant relationship stay if only OAs vs if only non-OAs
prevalent_ard_set = [code_chapter_to_field_mapper[x] for x in
                     list(one_firth_many_cox_mw_to_plot[one_firth_many_cox_mw_to_plot['Chapter 13'] == 1].index)]
ablation_test_dicts = []
for disease in prevalent_ard_set:
    ablation_test_dicts.append({
        'disease': disease,
        'only_oas': scipy.stats.mannwhitneyu(
            bma_similarities[(bma_similarities['disease1'] == disease) & (
                bma_similarities.disease2.isin(['p131868', 'p131870', 'p131872', 'p131874']))].bma_sem_sim,
            bma_similarities[(bma_similarities['disease1'] != disease) & (
                bma_similarities.disease2.isin(['p131868', 'p131870', 'p131872', 'p131874']))].bma_sem_sim,
            alternative='greater').pvalue,
        'only_non_oas': scipy.stats.mannwhitneyu(
            bma_similarities[(bma_similarities['disease1'] == disease) & (
                bma_similarities.disease2.isin(['p131900', 'p131964', 'p131972']))].bma_sem_sim,
            bma_similarities[(bma_similarities['disease1'] != disease) & (
                bma_similarities.disease2.isin(['p131900', 'p131964', 'p131972']))].bma_sem_sim,
            alternative='greater').pvalue
    })
ablation_test_df = pd.DataFrame(ablation_test_dicts)
ablation_test_long = ablation_test_df.melt(id_vars='disease', var_name='cox_ard_set', value_name='greater_pval')
ablation_test_long = significance_labelling(ablation_test_long, pvalue_col='greater_pval',
                                            fdr_permissive_group_by='cox_ard_set')

# plot -log10(raw pvals) of select ARD vs all chapter 13, only OAs, only non-OAs
ablation_test_to_plot = (one_firth_many_cox_mw[['firth_disease', 'chapter_13']]
                         .merge(ablation_test_df.rename(columns={'disease': 'firth_disease'}))
                         .melt(id_vars='firth_disease', var_name='cox_ard_set', value_name='greater_pval'))
ablation_test_to_plot['log_pval'] = -np.log10(ablation_test_to_plot['greater_pval'])
ablation_test_to_plot = ablation_test_to_plot.merge(
    disease_info[['code_chapter', 'disease_field']].rename(columns={'disease_field': 'firth_disease'})).drop(
    columns='firth_disease')
ablation_test_to_plot['cox_ard_set'] = ablation_test_to_plot['cox_ard_set'].replace(
    to_replace={'chapter_13': 'All Chapter 13',
                'only_oas': 'Only OAs',
                'only_non_oas': 'Only non-OAs'})
plt.figure(figsize=(8, 6))
sns.barplot(data=ablation_test_to_plot, x='code_chapter', y='log_pval', hue='cox_ard_set',
            hue_order=['All Chapter 13', 'Only OAs', 'Only non-OAs'])
plt.xlabel('Prevalent disease ICD-10 code (ICD-10 chapter)')
plt.gca().tick_params(axis='x', rotation=90)
plt.ylabel('-log10(p)')
plt.legend(title='Incident ARD set')
plt.axhline(y=-np.log10(0.05), color='grey', lw=0.8, ls='--', zorder=-1)
plt.tight_layout()
plt.savefig(
    f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/prevalent_ards_with_many_incident_c13_ards.png")
plt.savefig(
    f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/prevalent_ards_with_many_incident_c13_ards.svg")
plt.close()
