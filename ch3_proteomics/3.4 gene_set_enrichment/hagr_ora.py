"""
Bottom-up enrichment analyses of FDR-sig associated proteins in ageing/senescence genes,
 without first-degree interactors.
"""
import pandas as pd
from pyprind import ProgBar

from utils.significance_labelling import significance_labelling
from utils.enrichment_analyses import overrepresentation

import sys

""" Define task """
association = 'firth'

""" Read invariant data """
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Read in ontology data """
senescence_de_genes = pd.read_csv('~/data/external/hagr/signatures1.csv', sep=';')

ageing_de_genes = pd.read_excel('~/data/hagr/ageing_signatures.xlsx', sheet_name=None, header=1)

""" Read in protein data """
proteomics = pd.read_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{association}.csv")
sig_proteomics = proteomics[proteomics['fdr_sig'] == 1].copy()

""" Read in mapping data """
# read in protein -> UniProt conversions
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
senescence_overexp_genes = senescence_de_genes[senescence_de_genes.ovevrexp == 1].assign(annotation='senescence_overexp')
senescence_underexp_genes = senescence_de_genes[senescence_de_genes.underexp == 1].assign(annotation='senescence_underexp')

ageing_overexp_genes = ageing_de_genes['TableS3-Over.All'].assign(annotation='ageing_overexp').rename(columns={'Gene': 'gene_symbol'})
ageing_underexp_genes = ageing_de_genes['TableS7-Under.All'].assign(annotation='ageing_underexp').rename(columns={'Gene': 'gene_symbol'})

# restrict annotation background to background genes that have at least one relevant annotation
all_backgrounds = pd.concat([
    senescence_overexp_genes, senescence_underexp_genes, ageing_overexp_genes, ageing_underexp_genes
]).merge(all_qc_proteins[['gene_symbol', 'protein_field']])[['protein_field', 'annotation']]

""" Overrepresentation enrichment """
overrep_results = pd.DataFrame()
bar = ProgBar(len(sig_proteomics.disease.unique()), stream=sys.stdout, title='Overrepresentation enrichment')
for disease in sig_proteomics.disease.unique():
    protein_test_set = sig_proteomics[sig_proteomics.disease == disease].protein.to_list()

    overrep_result = overrepresentation(test_set=protein_test_set,
                                        background_annotations=all_backgrounds,
                                        entity_col_id='protein_field', annotation_col_id='annotation')
    overrep_results = pd.concat([overrep_results, overrep_result.assign(disease=disease)])

    bar.update()

# tidy
overrep_results = (overrep_results.reset_index(drop=True)
                   .merge(disease_info[['disease_field', 'icd10_three_letter']]
                          .rename(columns={'disease_field': 'disease'}))).drop_duplicates(ignore_index=True)

# significance labelling
overrep_results_labelled = significance_labelling(overrep_results, fdr_permissive_group_by='disease', alpha=0.05)
overrep_results_nom_sig = overrep_results_labelled[overrep_results_labelled['nom_sig'] == 1].copy()

# export
overrep_results.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/hagr_ora.csv",
    index=False)
overrep_results_nom_sig.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/hagr_ora_nom_sig.csv",
    index=False)

""" Disease-specific further analysis: Firth M81 """
m81_ageing_underexp = (ageing_underexp_genes
                       .merge(all_qc_proteins[['gene_symbol', 'protein_field']])
                       .merge(sig_proteomics[sig_proteomics.disease == 'p131964']
                              .rename(columns={'protein': 'protein_field'})))
