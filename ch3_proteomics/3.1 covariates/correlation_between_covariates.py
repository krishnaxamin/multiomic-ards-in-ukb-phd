"""
Examine correlations between the proteomics covariates.
"""
from pandas import read_csv, DataFrame, to_datetime
from collections import Counter
from pyprind import ProgBar
from scipy.stats import pearsonr
from scipy.stats.contingency import association
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys

matplotlib.use('TkAgg')

# import covariates
participant_technical_covars = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv').drop('p23099_i0', axis=1)
protein_covars = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv')
participant_lifestyle_covars = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv')

# cohort_eids = list(read_csv('ukbiobank/proteomics/proteomics_pan_ukbb_eur_eids.csv', sep='\s+')['eid'])

participant_covars = participant_lifestyle_covars.merge(participant_technical_covars, on='eid')

# create mean date sample processed
protein_covars['processing_start_datetime'] = to_datetime(protein_covars['Processing_StartDate'])
protein_covars_mean_dates = protein_covars.groupby('PlateID')['processing_start_datetime'].mean().reset_index().rename(columns={'PlateID': 'p30901_i0', 'processing_start_datetime': 'mean_processing_start_date'})

# merge into participant covariates and calculate storage time
covars = participant_covars.merge(protein_covars_mean_dates, on='p30901_i0')
covars['storage_time'] = (covars['mean_processing_start_date'] - to_datetime(covars['p53_i0'])).dt.days

# select covariates that would be used in analysis without pruning of highly collinear variables
# genetic PCs excluded for now - as need to be obtained via Apocrita (PCA file too big)
covars_for_analysis = (['p21003_i0', 'p31', 'p54_i0', 'storage_time', 'p1558_i0', 'p20116_i0', 'p23099_i0', 'p22189',
                       'p6138_i0', 'p884_i0', 'p1438_i0', 'p1458_i0', 'p1329_i0', 'grip_strength'] +
                       ['p1349_i0', 'p1359_i0', 'p1389_i0'] +
                       ['p' + str(x) + '_i0' for x in range(1309, 1329, 10)]) # + ['PC' + str(x) for x in range(1, 11)]
participant_covars_selected = covars[covars_for_analysis].copy()

""" Different tests for different combinations """
# continuous-continuous: Pearson's
# categorical-categorical: Cramer's V
# categorical-continuous -> dummy continuous-continuous: Pearson's

continuous_covars = ['p21003_i0', 'storage_time', 'p23099_i0', 'p22189', 'p884_i0', 'p1438_i0', 'p1458_i0', 'grip_strength', 'p1309_i0', 'p1319_i0']
covars_for_analysis_covar_type = {k: ('continuous' if k in continuous_covars else 'categorical') for k in covars_for_analysis}

covars_heatmap_v2 = DataFrame(columns=covars_for_analysis, index=covars_for_analysis)
covars_heatmap_v2_pval = DataFrame(columns=covars_for_analysis, index=covars_for_analysis)

covar_pairs = list(itertools.permutations(covars_for_analysis, 2))

bar = ProgBar(iterations=len(covar_pairs), stream=sys.stdout, title='Calculating correlations')
for pair in covar_pairs:
    x = pair[0]
    y = pair[1]
    x_type = covars_for_analysis_covar_type[x]
    y_type = covars_for_analysis_covar_type[y]

    if x_type == 'continuous' and y_type == 'continuous':
        # Pearson's
        pearson = pearsonr(participant_covars_selected[x], participant_covars_selected[y])
        covars_heatmap_v2.loc[x, y] = pearson.statistic
        covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
    elif x_type == 'categorical' and y_type == 'categorical':
        # Cramer's V
        covars_heatmap_v2.loc[x, y] = association(np.column_stack([
            participant_covars_selected[x].astype('category').cat.codes + 1,
            participant_covars_selected[y].astype('category').cat.codes + 1]))
    else:
        # Categorical to continuous, then Pearson's
        if x_type == 'categorical':
            pearson = pearsonr(participant_covars_selected[x].astype('category').cat.codes,
                               participant_covars_selected[y])
            covars_heatmap_v2.loc[x, y] = pearson.statistic
            covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
        elif y_type == 'categorical':
            pearson = pearsonr(participant_covars_selected[x],
                               participant_covars_selected[y].astype('category').cat.codes)
            covars_heatmap_v2.loc[x, y] = pearson.statistic
            covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
        else:
            print(f"Issue: {x} & {y}")
    bar.update()

# absolute the similarities because we want things that are anti-linear to still be considered as collinear
covars_heatmap_v2.to_csv('~/ch3_proteomics/3.1 covariates/pearsons_cramers_coeffs_between_covariates_heatmap.csv', index=False)
covars_heatmap_v2_abs_float = abs(covars_heatmap_v2.astype(np.float64))

# remove nominally insignificant Pearson's results
covars_heatmap_v2_abs_float_nom_sig = covars_heatmap_v2_abs_float[covars_heatmap_v2_pval.astype(np.float64).fillna(-1) < 0.05].fillna(0)
covars_heatmap_v2_nom_sig = covars_heatmap_v2[covars_heatmap_v2_pval.astype(np.float64).fillna(-1) < 0.05].fillna(0)

# precompute linkages based on the calculated similarity measures
linkages = linkage(squareform(1 - covars_heatmap_v2_abs_float_nom_sig.values, checks=False))

# clustermap using the linkages
g = sns.clustermap(covars_heatmap_v2_nom_sig, center=0, cmap='PiYG',
                   cbar_kws={'label': 'Pearson\'s rho/\nCramer\'s V'}, xticklabels=1, yticklabels=1,
                   row_linkage=linkages, col_linkage=linkages)
g.savefig('~/ch3_proteomics/3.1 covariates/plots/pearsons_cramers_coeffs_between_covariates_clustermap.svg')
g.savefig('~/ch3_proteomics/3.1 covariates/plots/pearsons_cramers_coeffs_between_covariates_clustermap.png')

# assess average similarity of the clusters - only include clusters whose mean similarity > 90th quantile (0.324)
np.quantile(covars_heatmap_v2_abs_float_nom_sig.values, 0.9)
# cluster 1: p31, grip_strength, p23099_i0 = 0.643
cluster1 = ['p31', 'grip_strength', 'p23099_i0']
abs(covars_heatmap_v2.loc[cluster1, cluster1].fillna(0)).values.sum() / 2 / ((len(cluster1) ** 2 - len(cluster1)) / 2)
# cluster 2: p6138_i0, p1329_i0, p1349_i0, p1359_i0, p54_i0, p20116_i0, p1558_i0, p1389_i0 = 0.318
cluster2 = ['p6138_i0', 'p1329_i0', 'p1349_i0', 'p1359_i0', 'p54_i0', 'p20116_i0', 'p1558_i0', 'p1389_i0']
abs(covars_heatmap_v2.loc[cluster2, cluster2].fillna(0)).values.sum() / 2 / ((len(cluster2) ** 2 - len(cluster2)) / 2)
# cluster 3: p1289_i0, p1299_i0, p1309_i0, p1319_i0, p1528_i0, p1458_i0 = 0.127
cluster3 = ['p1309_i0', 'p1319_i0', 'p1458_i0']
abs(covars_heatmap_v2.loc[cluster3, cluster3].fillna(0)).values.sum() / 2 / ((len(cluster3) ** 2 - len(cluster3)) / 2)

""" Variation of body fat % within sex """
plt.figure(figsize=(6, 8))
ax = sns.violinplot(data=participant_covars_selected, x='p31', y='p23099_i0', hue='p31')
plt.xlabel('Sex')
plt.ylabel('Body fat %')
plt.tight_layout()
plt.savefig('~/ch3_proteomics/3.1 covariates/plots/body_fat_perc_var_within_sex.png')
plt.savefig('~/ch3_proteomics/3.1 covariates/plots/body_fat_perc_var_within_sex.svg')
plt.close()
