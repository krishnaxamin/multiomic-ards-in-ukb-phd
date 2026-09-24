"""
Reactome GSEA of FDR-sig associated proteins ranked by log(n_disease) + 0.01 * median(beta).
Asking: are commonly associated proteins involved in any particular pathways?
"""
from utils.enrichment_analyses import slim_annotations

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
import numpy as np
import networkx as nx
import gseapy as gp

mpl.use('TkAgg')

""" Define task """
association = 'cox'

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
# if ontology == 'reactome':
reactome_annotations_full = pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None)
reactome_annotations_full = reactome_annotations_full[
    (reactome_annotations_full[6] == 'TAS') & (reactome_annotations_full[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                     axis=1).set_axis(
    ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()

""" Read in protein data """
data_sig = pd.read_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{association}.csv")
data_sig = data_sig[data_sig['fdr_sig'] == 1].copy()
data_sig['effect_direction'] = -1
data_sig.loc[data_sig['beta'] > 0, 'effect_direction'] = 1

""" Set up GSEA protein rankings """
# get number of diseases a protein is associated with and the median absolute beta for that protein across those diseases
protein_rankings = pd.DataFrame()
for protein in data_sig.protein.unique():
    data_protein = data_sig[data_sig['protein'] == protein].copy()
    protein_summaries = pd.concat([pd.DataFrame([{'protein': protein,
                                                  'count': len(df) * np.sign(effect_dir),
                                                  'median_beta': np.median(df.beta)}])
                                   for effect_dir, df, in data_protein.groupby('effect_direction')])
    protein_rankings = pd.concat([protein_rankings, protein_summaries])
protein_rankings.sort_values(by=['count', 'median_beta'], ignore_index=True, inplace=True)
protein_rankings['log_count'] = np.log(abs(protein_rankings['count'])) * np.sign(protein_rankings['count'])

# set a ranking metric based on the number of diseases and the median beta
protein_rankings['gsea_metric'] = protein_rankings['log_count'] + 0.01 * protein_rankings['median_beta']

""" Read in mapping data """
# read in protein -> UniProt conversions
olink_uniprot_mappings = pd.read_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv')
protein_rankings_uniprot = protein_rankings.merge(
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

""" Make gene sets """
# isolate relevant annotations
annotations_for_ranked_proteins = annotations_for_protein_background[
    annotations_for_protein_background['id'].isin(list(protein_rankings_uniprot['uniprot']))].copy()

# rudimentary check that slimming has worked
assert (-np.log(annotations_for_ranked_proteins.reactome_annotation_id.value_counts() / len(
    annotations_for_ranked_proteins.id.unique()))).max() < -np.log(10 / 19026)

# make gene set lists
gene_sets = {annotation_id: list(set(annotation_df['id'].to_list())) for annotation_id, annotation_df in
             annotations_for_ranked_proteins.groupby('reactome_annotation_id')}

""" Run GSEA """
rnk = protein_rankings_uniprot[['uniprot', 'gsea_metric']].copy()
rnk.index = rnk['uniprot']
rnk.drop(columns='uniprot', inplace=True)
rnk.columns = [0]
gsea_result = gp.prerank(rnk=rnk, gene_sets=gene_sets, threads=1, permutation_num=1000, seed=42, verbose=True,
                         min_size=1, max_size=10000, outdir=None)
gsea_result_ordered = gsea_result.res2d.sort_values('FDR q-val', ascending=True)
gsea_result_ordered_nompval = gsea_result.res2d.sort_values('NOM p-val', ascending=True)

""" Export results """
# add annotation names to results
gsea_results_to_export = gsea_result_ordered.merge(
    reactome_annotations_full[['reactome_annotation_id', 'annotation_name']].
    rename(columns={'reactome_annotation_id': 'Term'}).drop_duplicates(), on='Term')
gsea_results_to_export.to_csv(
    '~/data/internal/proteomics/enrichment_analysis/' + association + '/reactome_gsea.csv',
    index=False)

""" Plot results """
import math

gsea_results_to_export = pd.read_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_gsea.csv")
gsea_results_nom_sig = gsea_results_to_export[gsea_results_to_export['NOM p-val'] < 0.05].sort_values(by='NOM p-val')[
    ['annotation_name', 'NES']]

plt.figure(figsize=(10, 8))
sns.barplot(data=gsea_results_nom_sig, x='NES', y='annotation_name')
plt.gca().tick_params(axis='y', labelsize=8 - 1 * math.floor(len(gsea_results_nom_sig) / 60))
plt.axvline(x=0, color='black', lw=1, ls='-', zorder=3)
plt.grid(True, which='major', axis='x', color='lightgrey', lw=0.8, ls='--')
plt.gca().set_axisbelow(True)
plt.xlabel('Normalised enrichment score (NES)')
plt.ylabel('Reactome term')
plt.tight_layout()
plt.savefig(
    f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/common_proteins_gsea_reactome.png")
plt.savefig(
    f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/common_proteins_gsea_reactome.svg")
plt.close()

""" Parentage analysis """
from utils.enrichment_analyses import reactome_parentage_analysis

# positive NES
_ = reactome_parentage_analysis(annotation_set=gsea_results_to_export[
    (gsea_results_to_export['NOM p-val'] < 0.05) & (gsea_results_to_export.NES > 0)].Term.to_list(),
                                plot_path=f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/common_proteins_gsea_reactome_positiveNES_parentage_analysis",
                                export_plot=True)

# negative NES
_ = reactome_parentage_analysis(annotation_set=gsea_results_to_export[
    (gsea_results_to_export['NOM p-val'] < 0.05) & (gsea_results_to_export.NES < 0)].Term.to_list(),
                                plot_path=f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/common_proteins_gsea_reactome_negativeNES_parentage_analysis",
                                export_plot=True)
