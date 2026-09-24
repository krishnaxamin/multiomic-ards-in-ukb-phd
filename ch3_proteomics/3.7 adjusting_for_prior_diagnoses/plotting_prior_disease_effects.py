""" Plot the prior-disease effects calculated on phenotype data in R. """
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

pcs_type = 'before_ards'

""" Prior disease effects on phenotype """

# region Join together results from logistic regressions on Apocrita
qu10_50_disease_fields = list(pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])
results_list = []
for field in qu10_50_disease_fields:
    if pcs_type == 'before_blood_before_ards' and field == 'p131036':
        continue
    result = pd.read_csv(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/{pcs_type}_effect_on_phenotype/{field}_prior_disease_pc_effects_on_phenotype.csv")
    result['disease_field'] = field
    results_list.append(result)
results = pd.concat(results_list)
results.to_csv(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/prior_disease_pc_{pcs_type}_effects_on_phenotype.csv", index=False)
# endregion

results = pd.read_csv(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/prior_disease_pc_{pcs_type}_effects_on_phenotype.csv")
num_pcs = len(results['covar'].unique())  # get number of PCs
num_tests = len(results)  # preserve number of tests for Bonferroni correction
# results.drop('var_explained', axis=1, inplace=True)
results = results[abs(results['coefficient']) < 20].copy()  # remove wildly high coefficients
results['logp'] = results['pvalue'].apply(lambda x: -1 * np.log10(x))
results_significant = results[results['pvalue'] < (0.05 / num_tests)].copy()
results_significant_no_inf = results_significant[results_significant['logp'] != np.inf].copy()
results_significant_inf = results_significant[results_significant['logp'] == np.inf].copy()
results_insig = results[results['pvalue'] >= (0.05 / num_tests)].copy()

# make custom colour bar, with viridis centre, grey for values under limit, black for values over limit
colors = ['grey', *plt.cm.viridis.colors, 'orange']
colour_bar_range = plt.Normalize(results_significant_no_inf.logp.min(), results_significant_no_inf.logp.max())
cmap = mpl.colors.LinearSegmentedColormap.from_list('custom_cmap', colors, N=256)
colour_bar = mpl.cm.ScalarMappable(cmap=cmap, norm=colour_bar_range)

fig, ax = plt.subplots(figsize=[11.693, 8.268])
sns.stripplot(data=results_insig, x='covar', y='coefficient', color='grey', order=['PC' + str(x) for x in range(1, num_pcs + 1)])
sns.stripplot(data=results_significant_inf, x='covar', y='coefficient', color='orange',
              order=['PC' + str(x) for x in range(1, num_pcs + 1)])
sns.stripplot(data=results_significant_no_inf, x='covar', y='coefficient', hue='logp',
              order=['PC' + str(x) for x in range(1, num_pcs + 1)], palette='viridis')
ax.get_legend().remove()
ax.figure.colorbar(colour_bar, extend='both', ax=ax, label='-log10(p)')
ax.set_xlabel('Prior disease principal component', fontsize=15)
if len(results) < num_tests:
    ax.set_ylabel('Log-odds ratio (|coefficient| > 20 not displayed for ease of visual)', fontsize=15)
else:
    ax.set_ylabel('Log-odds ratio', fontsize=15)
ax.set_title('Log-odds ratios from a logistic regression fitting PCs from prior disease info \nto phenotype. Insignificant values are greyed out', fontsize=15)
plt.axhline(y=0, color='grey', linestyle='dashed')
plt.xticks(rotation=90)
fig.tight_layout()
fig.savefig(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/plots/prior_disease_pc_{pcs_type}_effects_on_phenotype.png")
fig.savefig(f"~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/plots/prior_disease_pc_{pcs_type}_effects_on_phenotype.svg")
plt.close()
