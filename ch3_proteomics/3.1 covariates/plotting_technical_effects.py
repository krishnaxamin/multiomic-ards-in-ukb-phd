""" Plot the covariate effects calculated on protein and phenotype data in R """
from pandas import read_csv, concat

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

""" Covariate effects on proteins """
# region Join together results from linear regressions on Apocrita
protein_fields = list(read_csv('~/data/internal/proteomics/protein_fields.csv')['protein_field'])
results_list = []
for field in protein_fields:
    result = read_csv('~/ch3_proteomics/3.1 covariates/technical_effects/on_proteins/' + field + '_technical_effects_on_proteins.csv')
    result['protein_field'] = field
    results_list.append(result)
results = concat(results_list)
results.to_csv('~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_proteins.csv', index=False)
# endregion

results = read_csv('~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_proteins.csv')
# results.drop('var_explained', axis=1, inplace=True)
results['var_explained_perc'] = results['var_explained'].apply(lambda x: 100 * x)
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['age', 'sex', 'centre', 'date_centre', 'plate', 'well', 'plate_row', 'plate_column', 'date_processing', 'storage_time', 'batch'],
                value=['Age', 'Sex', 'Centre', 'Date at \ncentre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \nprocessed', 'Storage \ntime', 'Batch'],
                inplace=True)
results_significant = results[results['pvalue'] < (0.05 / (len(results)))].copy()
results_significant_no_inf = results_significant[results_significant['logp'] != np.inf].copy()
results_significant_inf = results_significant[results_significant['logp'] == np.inf].copy()
results_insig = results[results['pvalue'] >= (0.05 / (len(results)))].copy()

# make custom colour bar, with viridis centre, grey for values under limit, black for values over limit
colors = ['grey', *plt.cm.viridis.colors, 'orange']
colour_bar_range = plt.Normalize(results_significant_no_inf.logp.min(), results_significant_no_inf.logp.max())
cmap = mpl.colors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)
colour_bar = mpl.cm.ScalarMappable(cmap=cmap, norm=colour_bar_range)

fig, ax = plt.subplots(figsize=[11.693, 8.268])

sns.stripplot(data=results_insig, x='covar', y='var_explained', color='grey',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_inf, x='covar', y='var_explained', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='var_explained', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \nprocessed', 'Storage \ntime', 'Batch'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('R-squared value', fontsize=15)
ax.set_title('R-squared values for models fitting covariates to protein data. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0.01, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_proteins.png')
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_proteins.svg')

""" Covariate effects on prior phenotype """

# region Join together results from logistic regressions on Apocrita
qu10_50_disease_fields = list(read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])
results_list = []
for field in qu10_50_disease_fields:
    if field == 'p131036':
        continue
    result = read_csv('~/ch3_proteomics/3.1 covariates/technical_effects/on_prior_phenotype/' + field + '_technical_effects_on_phenotype.csv')
    result['disease_field'] = field
    results_list.append(result)
results = concat(results_list)
results.to_csv('~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_prior_phenotype.csv', index=False)
# endregion

results = read_csv('~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_prior_phenotype.csv')
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['age', 'sex', 'date_centre', 'centre', 'plate', 'well', 'plate_row', 'plate_column', 'mean_date_processing', 'storage_time', 'batch'],
                value=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'],
                inplace=True)
results_significant = results[results['pvalue'] < (0.05 / (len(results)))].copy()
results_significant_no_inf = results_significant[results_significant['logp'] != np.inf].copy()
results_significant_inf = results_significant[results_significant['logp'] == np.inf].copy()
results_insig = results[results['pvalue'] >= (0.05 / (len(results)))].copy()

# make custom colour bar, with viridis centre, grey for values under limit, black for values over limit
colors = ['grey', *plt.cm.viridis.colors, 'orange']
colour_bar_range = plt.Normalize(results_significant_no_inf.logp.min(), results_significant_no_inf.logp.max())
cmap = mpl.colors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)
colour_bar = mpl.cm.ScalarMappable(cmap=cmap, norm=colour_bar_range)

fig, ax = plt.subplots(figsize=[11.693, 8.268])
sns.stripplot(data=results_insig, x='covar', y='coefficient', color='grey',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_inf, x='covar', y='coefficient', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='coefficient', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('Log-odds ratio', fontsize=15)
ax.set_title('Log-odds ratios from a logistic regression fitting covariates to prior phenotype. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_prior_phenotype.png')
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_prior_phenotype.svg')

""" Covariate effects on future phenotype """

results = read_csv('~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_future_phenotype.csv')
# results.drop('var_explained', axis=1, inplace=True)
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['age', 'sex', 'date_centre', 'centre', 'plate', 'well', 'plate_row', 'plate_column', 'mean_date_processing', 'storage_time', 'batch'],
                value=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'],
                inplace=True)
results_significant = results[results['pvalue'] < (0.05 / (len(results)))].copy()
results_significant_no_inf = results_significant[results_significant['logp'] != np.inf].copy()
results_significant_inf = results_significant[results_significant['logp'] == np.inf].copy()
results_insig = results[results['pvalue'] >= (0.05 / (len(results)))].copy()

# make custom colour bar, with viridis centre, grey for values under limit, black for values over limit
colors = ['grey', *plt.cm.viridis.colors, 'orange']
colour_bar_range = plt.Normalize(results_significant_no_inf.logp.min(), results_significant_no_inf.logp.max())
cmap = mpl.colors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)
colour_bar = mpl.cm.ScalarMappable(cmap=cmap, norm=colour_bar_range)

fig, ax = plt.subplots(figsize=[11.693, 8.268])
sns.stripplot(data=results_insig, x='covar', y='coefficient', color='grey',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_inf, x='covar', y='coefficient', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='coefficient', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Mean \ndate \nsample \nprocessed', 'Storage \ntime', 'Batch'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('Log-hazard ratio', fontsize=15)
ax.set_title('Log-hazard ratios from a Cox regression fitting covariates to future phenotype. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_future_phenotype.png')
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/technical_effects_on_future_phenotype.svg')
