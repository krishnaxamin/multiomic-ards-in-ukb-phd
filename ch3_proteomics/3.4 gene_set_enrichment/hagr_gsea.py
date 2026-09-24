"""
Testing commonly FDR-sig proteins for enrichment in ageing/senescence genes via GSEA.
Ranking metric = log(n_disease) + 0.01 * median(abs(beta)), stratified by effect direction.
"""
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
senescence_de_genes = pd.read_csv('~/data/external/hagr/signatures1.csv', sep=';')

ageing_de_genes = pd.read_excel('~/data/hagr/ageing_signatures.xlsx', sheet_name=None, header=1)

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
# read in protein <-> UniProt <-> gene name conversions
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
senescence_de_genes_olink = senescence_de_genes.merge(all_qc_proteins[['gene_symbol', 'protein_field']])
senescence_overexp_genes = senescence_de_genes_olink[senescence_de_genes_olink['ovevrexp'] == 1].copy()
senescence_underexp_genes = senescence_de_genes_olink[senescence_de_genes_olink['underexp'] == 1].copy()

ageing_overexp_genes = ageing_de_genes['TableS3-Over.All'].merge(
    all_qc_proteins[['gene_symbol', 'protein_field']].rename(columns={'gene_symbol': 'Gene'}))
ageing_underexp_genes = ageing_de_genes['TableS7-Under.All'].merge(
    all_qc_proteins[['gene_symbol', 'protein_field']].rename(columns={'gene_symbol': 'Gene'}))

""" Make gene sets """
gene_sets = {'senescence_overexp': senescence_overexp_genes.protein_field.to_list(),
             'senescence_underexp': senescence_underexp_genes.protein_field.to_list(),
             'ageing_overexp': ageing_overexp_genes.protein_field.to_list(),
             'ageing_underexp': ageing_underexp_genes.protein_field.to_list()}

""" Run GSEA """
rnk = protein_rankings[['protein', 'gsea_metric']].copy()
rnk.index = rnk['protein']
rnk.drop(columns='protein', inplace=True)
rnk.columns = [0]
gsea_result = gp.prerank(rnk=rnk, gene_sets=gene_sets, threads=1, permutation_num=1000, seed=42, verbose=True,
                         min_size=1, max_size=10000, outdir=None)
gsea_result_ordered = gsea_result.res2d.sort_values('FDR q-val', ascending=True)
gsea_result_ordered_nompval = gsea_result.res2d.sort_values('NOM p-val', ascending=True)

""" Export results """
gsea_result_ordered.to_csv(
    '~/data/internal/proteomics/enrichment_analysis/' + association + '/hagr_gsea.csv',
    index=False)

""" Run GSEA on HAGR without differentiating DE direction """
gene_sets_no_dir = {'senescence': [x.lower() for x in senescence_de_genes.gene_symbol.to_list()],
                    'ageing': [x.lower() for x in pd.concat(
                        [ageing_de_genes['TableS3-Over.All'], ageing_de_genes['TableS7-Under.All']]).Gene.to_list()]}

rnk = protein_rankings[['protein', 'gsea_metric']].copy()
rnk.index = rnk['protein']
rnk.drop(columns='protein', inplace=True)
rnk.columns = [0]
gsea_result_no_dir = gp.prerank(rnk=rnk, gene_sets=gene_sets_no_dir, threads=1, permutation_num=1000, seed=42,
                                verbose=True,
                                min_size=1, max_size=10000, outdir=None)
gsea_result_no_dir_ordered = gsea_result_no_dir.res2d.sort_values('FDR q-val', ascending=True)
gsea_result_no_dir_ordered_nompval = gsea_result_no_dir.res2d.sort_values('NOM p-val', ascending=True)
gsea_result_no_dir_ordered.to_csv(
    '~/data/internal/proteomics/enrichment_analysis/' + association + '/hagr_gsea_no_dir.csv',
    index=False)

""" Plot results (Firth and Cox separately) """
hagr_to_plot = pd.concat(
    [pd.read_csv(f"~/data/internal/proteomics/enrichment_analysis/{association}/hagr_gsea.csv"),
     pd.read_csv(
         f"~/data/internal/proteomics/enrichment_analysis/{association}/hagr_gsea_no_dir.csv")])
hagr_to_plot = (hagr_to_plot[['Term', 'NES', 'FDR q-val']]
                .replace(to_replace='senescence_overexp', value='Senescence (+)')
                .replace(to_replace='senescence_underexp', value='Senescence (-)')
                .replace(to_replace='ageing_overexp', value='Ageing (+)')
                .replace(to_replace='ageing_underexp', value='Ageing (-)')
                .replace(to_replace='senescence', value='Senescence (ALL)')
                .replace(to_replace='ageing', value='Ageing (ALL)')
                .rename(columns={'FDR q-val': 'FDR < 0.05', 'Term': 'Gene set'}))
hagr_to_plot['FDR < 0.05'] = (hagr_to_plot['FDR < 0.05'] < 0.05).astype(str)

plt.figure(figsize=(8, 6))
sns.barplot(data=hagr_to_plot, x='NES', y='Gene set', hue='FDR < 0.05', hue_order=['True', 'False'],
            palette=(plt.cm.coolwarm(1.0), plt.cm.coolwarm(0.0)))
plt.xlabel('Normalised enrichment score (NES)')
plt.grid(True, which='major', axis='x', color='lightgrey', lw=0.8, ls='--')
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig(f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/hagr_gsea.png")
plt.savefig(f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/hagr_gsea.svg")
plt.close()
