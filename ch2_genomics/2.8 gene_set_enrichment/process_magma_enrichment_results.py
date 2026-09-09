""" Process MAGMA enrichment results on SAIGE results. """
from statsmodels.stats.multitest import fdrcorrection
from itertools import combinations
from scipy.stats import fisher_exact
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from utils.significance_labelling import significance_labelling
from utils.mainali_alpha import mainali_alpha
from utils.enrichment_analyses import slim_annotations, propagation_up_ontology_hierarchy, overrepresentation

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

""" Set up """
adjustment = 'lifestyles'
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Gene set analysis - Reactome """
# get results files
magma_results_files = [x.name for x in
                       os.scandir('~/data/internal/genomics/magma/' + adjustment + '/enrichment_results/') if
                       x.is_file() and '.gsa.out' in x.name]
results = pd.DataFrame()
for results_file in magma_results_files:
    result = pd.read_csv('~/data/internal/genomics/magma/' + adjustment + '/enrichment_results/' + results_file,
                      sep='\s+', comment='#')
    disease = results_file.split('_')[0]
    if 'lifestyles' in disease:
        disease = disease.split('_')[0]
    result['disease'] = disease
    results = pd.concat([results, result])
results.reset_index(drop=True, inplace=True)
# add names - data accessed: 2025-06-10
reactome = pd.read_csv('~/data/external/reactome/NCBI2Reactome_PE_All_Levels.txt',
                    sep='\t', header=None)
reactome.columns = ['entrez_id', 'reactome_id', 'reactome_name', 'reactome_pathway_id', 'reactome_pathway_link',
                    'reactome_pathway_name', 'evidence', 'species']
reactome_human_tas = reactome[(reactome['species'] == 'Homo sapiens') & (reactome['evidence'] == 'TAS')].copy()
reactome_human_tas_for_merging = reactome_human_tas[['reactome_pathway_id', 'reactome_pathway_name']].rename(
    columns={'reactome_pathway_id': 'VARIABLE'}).drop_duplicates(ignore_index=True)
results_labelled = pd.merge(results, reactome_human_tas_for_merging, on='VARIABLE',
                         how='left').drop_duplicates(ignore_index=True)
results_labelled.to_csv('~/data/internal/genomics/magma/' + adjustment + '/reactome_enrichment_results_full.csv', index=False)

# annotate with significance labels
# permissive FDR and Bonferroni (correcting within disease only)
results_assoc_labelled = significance_labelling(results_labelled, pvalue_col='P', fdr_permissive_group_by='disease')
results_assoc_labelled = results_assoc_labelled[results_assoc_labelled['nom_sig'] == 1].copy()
results_assoc_labelled.to_csv(
    '~/data/internal/genomics/magma/' + adjustment + '/reactome_enrichment_results_assoc.csv', index=False)

# parentage analysis
from utils.enrichment_analyses import reactome_parentage_analysis
# E83
_ = reactome_parentage_analysis(annotation_set=results_assoc_labelled[(results_assoc_labelled.fdr_sig == 1) &
                                                                      (results_assoc_labelled.disease == 'p130820')].VARIABLE.to_list(),
                                plot_path=f"~/ch2_genomics/2.8 gene_set_enrichment/plots/{adjustment}/parentage_analysis_E83",
                                export_plot=True)
# I25
_ = reactome_parentage_analysis(annotation_set=results_assoc_labelled[(results_assoc_labelled.fdr_sig == 1) &
                                                                      (results_assoc_labelled.disease == 'p131306')].VARIABLE.to_list(),
                                plot_path=f"~/ch2_genomics/2.8 gene_set_enrichment/plots/{adjustment}/parentage_analysis_I25",
                                export_plot=True)

# all assoc terms
_ = reactome_parentage_analysis(annotation_set=list(results_assoc_labelled[results_assoc_labelled.fdr_sig == 1].VARIABLE.unique()),
                                plot_path=f"~/ch2_genomics/2.8 gene_set_enrichment/plots/{adjustment}/parentage_analysis_all_sig_terms",
                                export_plot=True)