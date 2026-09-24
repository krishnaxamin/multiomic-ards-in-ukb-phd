""" Assess disease similarity using Resnik BMA semantic similarity scores between top 20 (by effect) annotations. """
import itertools
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
import re
import scipy

from pyprind import ProgBar
from utils.enrichment_analyses import prep_annotation_files
from utils.semantic_similarity import resnik_bma_between_two_annotation_sets

mpl.use('TkAgg')

""" Define task """
association = 'cox'

""" Read invariant data """

# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
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
data = pd.read_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora.csv")
# impute infinite OR to be the maximum finite OR
data.loc[data['odds_ratio'] == np.inf, 'odds_ratio'] = data[data['odds_ratio'] != np.inf]['odds_ratio'].max()

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
#  (1) IC is dependent on genes, not terms, in the data. Slimming affects terms only, not genes
#  (2) slimming would remove high-level MICAs necessary for SemSim calcs
# annotations_for_protein_background = slim_annotations(annotations=annotations_for_protein_background,
#                                                       annotation_col_id='annotation_id', entity_col_id='id')

# calculate IC for all annotations FOR ONLY THOSE PROTEINS ANALYSED
annotations_ic = -np.log(annotations_for_protein_background.annotation_id.value_counts() / len(
    annotations_for_protein_background.id.unique()))

# no need to check that slimming has worked, since the background annotation set hasn't been slimmed
# if wanting to check that slimming is OK, filter to the pathways present in the enrichment results
# annotations_ic[data.annotation].max() gives the max possible SemSim for this dataset, obtained if every pathway in
#  both datasets is the most informative post-slimming pathway annotated to the assoc proteins
assert annotations_ic[data.annotation].max() < -np.log(10 / 19026)

""" Semantic-based similarity """
# calculate weights for all annotation-disease pairs
data['weights'] = -np.log10(data['pval']) * np.log10(data['odds_ratio'])

# get absolute value of log10(OR) - used for ranking the top most involved pathways
data['abs_log_or'] = abs(np.log10(data['odds_ratio']))

# for each disease group, for each combination of annotations across the diseases, calculate a similarity metric
disease_pairs = list(itertools.combinations(sorted(list(data['disease'].unique())), 2))
disease_annotation_sets = {disease: d.sort_values('abs_log_or', ascending=False, ignore_index=True).loc[:19,
                                    ['annotation']]['annotation'].to_list() for disease, d in data.groupby('disease')}

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

""" Label the MICAs with annotation names """
reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                header=None, names=['id', 'name', 'species'])
reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
bma_mica_contributions_export = bma_mica_contributions.merge(
    reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

""" Export semantic similarities """
bma_similarities.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_disease_resnik_bma.csv",
    index=False)
bma_mica_contributions_export.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_disease_resnik_micas.csv",
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


semsims_within_between_chapters_mwu(bma_similarities)

""" Firth-Firth and Cox-Cox heatmap """


def plot_intra_algo_heatmaps():

    firthfirth = pd.read_csv(
        f"~/data/internal/proteomics/pathway_based_similarities/firth_disease_resnik_bma.csv")
    coxcox = pd.read_csv(f"~/data/internal/proteomics/pathway_based_similarities/cox_disease_resnik_bma.csv")

    data_to_plot = pd.concat([firthfirth[['disease1', 'disease2', 'bma_sem_sim']],
                              coxcox[['disease1', 'disease2', 'bma_sem_sim']]
                             .rename(columns={'disease1': 'disease2', 'disease2': 'disease1'})])
    data_to_plot = (data_to_plot
                    .merge(disease_info[['disease_field', 'code_chapter']]
                           .rename(columns={'disease_field': 'disease1'}))
                    .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                    .merge(disease_info[['disease_field', 'code_chapter']]
                           .rename(columns={'disease_field': 'disease2'}))
                    .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))

    data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                              categories=sorted(data_to_plot['disease1'].unique()))
    data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                              categories=sorted(data_to_plot['disease2'].unique()))

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

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.scatterplot(data_to_plot, x='disease1_as_num', y='disease2_as_num', hue='bma_sem_sim',
                    palette='viridis', zorder=3, ax=ax, s=50)
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', title='Resnik BMA', fontsize=8)
    plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=9)
    plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=9)
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
    norm = mpl.colors.Normalize(vmin=data_to_plot.bma_sem_sim.min(),
                                vmax=data_to_plot.bma_sem_sim.max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # required boilerplate

    # remove legend, attach colorbar
    ax.get_legend().remove()
    plt.colorbar(sm, ax=ax, label='Semantic similarity')

    plt.tight_layout()
    plt.savefig(
        f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firthfirth_coxcox_heatmap_disease_pairs_by_resnik_bma.svg")
    plt.savefig(
        f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/firthfirth_coxcox_heatmap_disease_pairs_by_resnik_bma.png")
    plt.close()


plot_intra_algo_heatmaps()
