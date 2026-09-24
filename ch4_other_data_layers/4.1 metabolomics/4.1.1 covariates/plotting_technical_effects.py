""" Plot the covariate effects calculated on metabolite and phenotype data in R """
from pandas import read_csv, concat

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

""" Covariate effects on metabolites """
# region Join together results from linear regressions on Apocrita
metabolite_fields = list(read_csv('~/data/internal/metabolomics/metabolite_fields.csv')['metabolite_field'])
results_list = []
for field in metabolite_fields:
    result = read_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/on_metabolites/' + field + '_technical_effects_on_metabolites.csv')
    result['metabolite_field'] = field
    results_list.append(result)
results = concat(results_list)
results.to_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/technical_effects_on_metabolites.csv', index=False)
# endregion

results = read_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/technical_effects_on_metabolites.csv')
results['var_explained_perc'] = results['var_explained'].apply(lambda x: 100 * x)
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['Date', 'Storage time', 'Prepped-for time', 'Plate row', 'Plate column', 'Sample Measured Date', 'Sample Prepared Date'],
                            value=['Date at \ncentre', 'Storage \ntime', 'Prepped-for \ntime', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \n measured', 'Date \nsample \n prepared'],
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
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'])
sns.stripplot(data=results_significant_inf, x='covar', y='var_explained', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='var_explained', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('R-squared value', fontsize=15)
ax.set_title('R-squared values for models fitting covariates to metabolite data. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0.01, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/technical_effects_on_metabolites.png')
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/technical_effects_on_metabolites.svg')

""" Covariate effects on prior phenotype """

# region Join together results from logistic regressions on Apocrita
qu10_50_disease_fields = list(read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])
results_list = []
for field in qu10_50_disease_fields:
    result = read_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/on_prior_phenotype/' + field + '_technical_effects_on_phenotype.csv')
    result['disease_field'] = field
    results_list.append(result)
results = concat(results_list)
results.to_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/technical_effects_on_prior_phenotype.csv', index=False)
# endregion

results = read_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/technical_effects/technical_effects_on_prior_phenotype.csv')
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results.replace(to_replace=['Date', 'Storage time', 'Prepped-for time', 'Plate row', 'Plate column', 'Sample Measured Date', 'Sample Prepared Date'],
                            value=['Date at \ncentre', 'Storage \ntime', 'Prepped-for \ntime', 'Plate \nrow', 'Plate \ncolumn', 'Date \nsample \n measured', 'Date \nsample \n prepared'],
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
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'])
sns.stripplot(data=results_significant_inf, x='covar', y='coefficient', color='orange',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'])
sns.stripplot(data=results_significant_no_inf, x='covar', y='coefficient', hue='logp',
              order=['Age', 'Sex', 'Date at \ncentre', 'Centre', 'Plate', 'Spectrometer', 'Well', 'Plate \nrow', 'Plate \ncolumn', 'Prepped-for \ntime', 'Storage \ntime', 'Date \nsample \n prepared', 'Date \nsample \n measured'],
              palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Covariate', fontsize=15)
ax.set_ylabel('Log-odds ratio', fontsize=15)
ax.set_title('Log-odds ratios from a logistic regression fitting covariates to prior phenotype. \nInsignificant values are greyed out', fontsize=15)
plt.axhline(y=0, color='grey', linestyle='dashed')
fig.tight_layout()
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/technical_effects_on_prior_phenotype.png')
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/technical_effects_on_prior_phenotype.svg')
