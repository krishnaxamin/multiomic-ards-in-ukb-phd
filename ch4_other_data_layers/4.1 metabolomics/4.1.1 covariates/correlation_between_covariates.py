"""
Examine correlations between the metabolomics lifestyles and covariates.
"""
from pandas import read_csv, DataFrame, to_datetime, get_dummies
from collections import Counter
from pyprind import ProgBar
from scipy.stats import pearsonr, chi2_contingency
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


def entropy(labels):
    """Calculate the entropy of a list of labels."""
    n_labels = len(labels)

    if n_labels <= 1:
        return 0

    counts = np.array(list(Counter(labels).values()))  # counts of each unique value
    probs = counts / n_labels
    return -np.sum(probs * np.log(probs))  # entropy calculation: -1 * sum(probs * log(probs))


def conditional_entropy(x, y):
    """Calculate the conditional entropy of y given x."""
    labels_x, counts_x = np.unique(x, return_counts=True)
    n_samples = len(x)
    cond_entropy = 0.0

    for label_x, count_x in zip(labels_x, counts_x):
        indices = np.where(x == label_x)[0]  # get where X = x
        cond_entropy += (count_x / n_samples) * entropy(y[indices])  # p(x) * H(Y|(X=x)) [entropy of Y where X = x]. This representation of conditional entropy is valid for discrete random variables

    return cond_entropy


def theils_u(x, y):
    """Calculate Theil's U coefficient for variables x and y: y given x."""
    entropy_y = entropy(y)  # H(Y)
    cond_entropy_xy = conditional_entropy(x, y)  # H(Y|X)
    if entropy_y == 0:
        return 1
    else:
        u = (entropy_y - cond_entropy_xy) / entropy_y
        return u


# import covariates
technical_covars = read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_covars.csv').drop('p23099_i0', axis=1)
lifestyle_covars = read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_lifestyles.csv')

covars = lifestyle_covars.merge(technical_covars, on='eid')

# select covariates that would be used in analysis without pruning of highly collinear variables
# genetic PCs excluded for now - as need to be obtained via Apocrita (PCA file too big)
covars_for_analysis = (['p21003_i0', 'p31', 'p54_i0', 'storage_time', 'p1558_i0', 'p20116_i0', 'p23099_i0', 'p1160_i0',
                       'p22189', 'p884_i0', 'p1438_i0', 'p6138_i0', 'p1528_i0', 'p1458_i0', 'grip_strength'] +
                       ['p' + str(x) + '_i0' for x in range(1289, 1399, 10)]) # + ['PC' + str(x) for x in range(1, 11)]
covars_selected = covars[covars_for_analysis].copy()

""" Different tests for different combinations """
# continuous-continuous: Pearson's
# categorical-categorical: Cramer's V
# categorical-continuous -> dummy continuous-continuous: Pearson's

continuous_covars = ['p21003_i0', 'storage_time', 'p23099_i0', 'p22189', 'p884_i0', 'p1438_i0', 'p1458_i0', 'p1528_i0', 'grip_strength', 'p1289_i0', 'p1299_i0', 'p1309_i0', 'p1319_i0']
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
        pearson = pearsonr(covars_selected[x], covars_selected[y])
        covars_heatmap_v2.loc[x, y] = pearson.statistic
        covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
    elif x_type == 'categorical' and y_type == 'categorical':
        # Cramer's V
        covars_heatmap_v2.loc[x, y] = association(np.column_stack([
            covars_selected[x].astype('category').cat.codes + 1,
            covars_selected[y].astype('category').cat.codes + 1]))
    else:
        # Categorical to continuous, then Pearson's
        if x_type == 'categorical':
            pearson = pearsonr(covars_selected[x].astype('category').cat.codes,
                               covars_selected[y])
            covars_heatmap_v2.loc[x, y] = pearson.statistic
            covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
        elif y_type == 'categorical':
            pearson = pearsonr(covars_selected[x],
                               covars_selected[y].astype('category').cat.codes)
            covars_heatmap_v2.loc[x, y] = pearson.statistic
            covars_heatmap_v2_pval.loc[x, y] = pearson.pvalue
        else:
            print(f"Issue: {x} & {y}")
    bar.update()

covars_heatmap_v2.to_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/pearsons_cramers_coeffs_between_covariates_heatmap.csv', index=False)
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
g.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/pearsons_cramers_coeffs_between_covariates_clustermap.svg')
g.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/pearsons_cramers_coeffs_between_covariates_clustermap.png')

# cluster 1: p31, grip_strength, p23099_i0
cluster1 = ['p31', 'grip_strength', 'p23099_i0']
abs(covars_heatmap_v2.loc[cluster1, cluster1].fillna(0)).values.sum() / 2 / ((len(cluster1) ** 2 - len(cluster1)) / 2)
# cluster 2: p6138_i0, p1329_i0, p1349_i0, p1359_i0, p1369_i0, p54_i0, p20116_i0, p1558_i0, p1160_i0, p1379_i0, p1389_i0
cluster2 = ['p6138_i0', 'p1329_i0', 'p1349_i0', 'p1359_i0', 'p1369_i0', 'p54_i0', 'p20116_i0', 'p1558_i0', 'p1160_i0', 'p1379_i0', 'p1389_i0']
abs(covars_heatmap_v2.loc[cluster2, cluster2].fillna(0)).values.sum() / 2 / ((len(cluster2) ** 2 - len(cluster2)) / 2)
# cluster 3: p1289_i0, p1299_i0, p1309_i0, p1319_i0, p1528_i0, p1458_i0
cluster3 = ['p1289_i0', 'p1299_i0', 'p1309_i0', 'p1319_i0', 'p1528_i0', 'p1458_i0']
abs(covars_heatmap_v2.loc[cluster3, cluster3].fillna(0)).values.sum() / 2 / ((len(cluster3) ** 2 - len(cluster3)) / 2)
