""" Cluster ARDs using pathway-based similarities, after performing prerequisite Reactome ORA and semantic similarity calculations. """
import itertools
import pandas as pd
import sys
import numpy as np
import scipy
import statsmodels.stats.multitest as smm
import networkx as nx

from pyprind import ProgBar
from utils.enrichment_analyses import prep_annotation_files, overrepresentation, slim_annotations
from utils.semantic_similarity import resnik_bma_between_two_annotation_sets
from utils.significance_labelling import significance_labelling


""" Define task """
association = 'cox'

""" Read invariant data """

# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
field_to_code_mapper = dict(zip(disease_info.disease_field, disease_info.icd10_three_letter))
code_chapter_to_field_mapper = dict(zip(disease_info.code_chapter, disease_info.disease_field))

""" Read in ontology data """
reactome_annotations_full = pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None)
reactome_annotations_full = reactome_annotations_full[
    (reactome_annotations_full[6] == 'TAS') & (reactome_annotations_full[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                     axis=1).set_axis(
    ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()

""" Read in protein data """
proteomics = pd.read_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{association}_prior_disease.csv")
sig_proteomics = proteomics[proteomics['fdr_sig'] == 1].copy()

""" Read in mapping data """
# read in protein -> UniProt conversions
olink_uniprot_mappings = pd.read_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv')
sig_proteomics_uniprot = sig_proteomics.merge(
    olink_uniprot_mappings[['uniprot', 'protein_field']].rename(columns={'protein_field': 'protein'}))

""" Define protein and annotation background """
# --- protein background ---
# use the mapping file, as this contains all the UKB proteins before QC
# remove the 3 proteins calculated as high missingness during QC
high_missingness_proteins = pd.read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/proteomics/high_missingness_protein_fields.txt')
all_qc_proteins = olink_uniprot_mappings[
    ~olink_uniprot_mappings['protein_field'].isin(high_missingness_proteins['protein_field'])].copy()
protein_background = all_qc_proteins['uniprot'].to_list()

# --- annotation background ---
annotations_for_protein_background = reactome_annotations_full[
    reactome_annotations_full['id'].isin(protein_background)].copy()
annotations_for_protein_background = annotations_for_protein_background[
    ['id', 'reactome_annotation_id', 'annotation_name']].drop_duplicates()

# slim
annotations_for_protein_background = slim_annotations(annotations=annotations_for_protein_background,
                                                      annotation_col_id='reactome_annotation_id', entity_col_id='id')

""" Annotate all assoc proteins """
sig_proteomics_with_annotations = annotations_for_protein_background[
    annotations_for_protein_background['id'].isin(sig_proteomics_uniprot.uniprot.to_list())].copy()

# rudimentary check that slimming has worked
assert (-np.log(sig_proteomics_with_annotations.reactome_annotation_id.value_counts() / len(
    sig_proteomics_with_annotations.id.unique()))).max() < -np.log(10 / 19026)

""" Overrepresentation enrichment """
overrep_results = pd.DataFrame()
bar = ProgBar(len(sig_proteomics.disease.unique()), stream=sys.stdout, title='Overrepresentation enrichment')
for disease in sig_proteomics.disease.unique():
    protein_test_set = sig_proteomics_uniprot[sig_proteomics_uniprot.disease == disease].uniprot.to_list()

    overrep_result = overrepresentation(test_set=protein_test_set,
                                        background_annotations=annotations_for_protein_background,
                                        entity_col_id='id', annotation_col_id='reactome_annotation_id')
    overrep_results = pd.concat([overrep_results, overrep_result.assign(disease=disease)])

    bar.update()

# tidy
overrep_results = (overrep_results.reset_index(drop=True)
                   .merge(reactome_annotations_full[['reactome_annotation_id', 'annotation_name']]
                          .rename(columns={'reactome_annotation_id': 'annotation'}))
                   .merge(disease_info[['disease_field', 'icd10_three_letter']]
                          .rename(columns={'disease_field': 'disease'}))).drop_duplicates(ignore_index=True)

# significance labelling
overrep_results_labelled = significance_labelling(overrep_results, fdr_permissive_group_by='disease', alpha=0.05)
overrep_results_nom_sig = overrep_results_labelled[overrep_results_labelled['nom_sig'] == 1].copy()

# export
overrep_results.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora_prior_disease.csv",
    index=False)
overrep_results_nom_sig.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora_prior_disease_nom_sig.csv",
    index=False)

""" Semantic similarity"""
""" Annotations for SemSim calculations """
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
hierarchy = pd.read_csv(
    '~/data/external/reactome/reactome_id_ancestor_relations.csv').set_axis(
    ['id', 'ancestors'], axis=1)

""" Read in data to use to assess commonalities """
data = pd.read_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora_prior_disease.csv")
# impute infinite OR to be the maximum finite OR
data.loc[data['odds_ratio'] == np.inf, 'odds_ratio'] = data[data['odds_ratio'] != np.inf]['odds_ratio'].max()

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

""" Calculate semantic similarity """
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

""" Export semantic similarities """
bma_similarities.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_disease_resnik_bma_prior_disease.csv",
    index=False)

""" Louvain clustering """

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
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_resnik_bma_disease_groups_prior_disease.csv",
    index=False)
