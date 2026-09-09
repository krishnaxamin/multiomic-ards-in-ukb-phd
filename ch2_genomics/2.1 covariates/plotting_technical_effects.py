""" Plot results from technical covariate effects on phenotype using logistic regression """
from pandas import read_csv, concat

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

# region Join together results from logistic regressions on Apocrita
qu10_50_disease_fields = list(read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])
results_list = []
for field in qu10_50_disease_fields:
    result = read_csv('~/ch2_genomics/2.1 covariates/technical_effects/' + field + '_technical_effects_on_phenotype.csv')
    result['disease_field'] = field
    results_list.append(result)
results = concat(results_list)
results.to_csv('~/ch2_genomics/2.1 covariates/technical_effects_on_phenotype.csv', index=False)
# endregion

results = read_csv('~/ch2_genomics/2.1 covariates/technical_effects_on_phenotype.csv')
# results.drop('var_explained', axis=1, inplace=True)
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['Date', 'Plate row', 'Plate column'],
                value=['Date at \ncentre', 'Plate \nrow', 'Plate \ncolumn'],
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
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Batch'])
sns.stripplot(data=results_significant_inf, x='covar', y='coefficient', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Batch'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='coefficient', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Batch'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('Log-odds ratio', fontsize=15)
ax.set_title('Log-odds ratios from a logistic regression fitting covariates to phenotype. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch2_genomics/2.1 covariates/plots/technical_effects_on_phenotype.png')
fig.savefig('~/ch2_genomics/2.1 covariates/plots/technical_effects_on_phenotype.svg')
