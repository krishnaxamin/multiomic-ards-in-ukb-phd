""" Process MAGMA gene-level results on SAIGE results. """
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

""" Gene analysis """
# get results files
magma_results_files = [x.name for x in os.scandir('~/data/internal/genomics/magma/' + adjustment + '/gene_level_results/') if
                       x.is_file() and '.genes.out' in x.name]
results = pd.DataFrame()
for results_file in magma_results_files:
    result = pd.read_csv('~/data/internal/genomics/magma/' + adjustment + '/gene_level_results/' + results_file, sep='\s+')
    disease = results_file.split('.')[0]
    if 'lifestyles' in disease:
        disease = disease.split('_')[0]
    result['disease'] = disease
    results = pd.concat([results, result])
results.reset_index(drop=True, inplace=True)
results.to_csv('~/data/internal/genomics/magma/' + adjustment + '/gene_level_results_full.csv', index=False)

# extract invariant gene info
gene_info = results.drop(['NSNPS', 'NPARAM', 'ZSTAT', 'P', 'disease'], axis=1).drop_duplicates()
# len(gene_info.GENE.unique())
gene_info.to_csv('~/data/internal/genomics/magma/tested_genes.csv', index=False)

# get slimmed MAGMA results
results_slimmed = results[['GENE', 'ZSTAT', 'P', 'disease']].copy()
# annotate with significance labels
results_slimmed.loc[results_slimmed['P'] < 0.05, 'nom_sig'] = 1
results_slimmed['fdr_sig'], results_slimmed['pval_fdr_corrected'] = fdrcorrection(results_slimmed['P'], alpha=0.05)
results_slimmed.loc[results_slimmed['P'] < 0.05 / len(results_slimmed), 'bonf_sig'] = 1
# permissive FDR (correcting within disease only)
results_slimmed_assoc_labelled = significance_labelling(results_slimmed, pvalue_col='P',
                                                        fdr_permissive_group_by='disease')
results_slimmed_assoc_labelled = results_slimmed_assoc_labelled[results_slimmed_assoc_labelled['nom_sig'] == 1].copy()
results_slimmed_assoc_labelled.to_csv('~/data/internal/genomics/magma/' + adjustment + '/gene_level_results_assoc.csv',
                                      index=False)
