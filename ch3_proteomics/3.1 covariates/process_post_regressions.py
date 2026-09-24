"""
Collate post-regression protein residuals and assess the effect each regression has on variance explained by the lifestyle factors.
"""
from pandas import read_csv, merge, concat, DataFrame
from datetime import datetime
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, robust_scale
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from pyprind import ProgBar

import sys
import matplotlib.pyplot as plt
import seaborn as sns

input_data = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_data.csv')
protein_fields = list(input_data.columns)[1:]

# merging residuals
residuals = input_data[['eid']].copy()
bar = ProgBar(len(protein_fields), stream=sys.stdout, title='Collating per-protein final residuals')
for protein in protein_fields:
    residual = read_csv('~/ch3_proteomics/3.1 covariates/regressing_out/post_regression_protein_data/' + protein + '.csv')
    residuals = merge(residuals, residual, how='left')
    bar.update()

# robust_scale the post-regressions data and export
scaled_residuals = residuals.copy()
scaled_residuals.iloc[:, 1:] = robust_scale(residuals.iloc[:, 1:])
scaled_residuals.to_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_post_regressions_robust_scaled_data.csv', index=False)

# vars explained at each step
vars_explained_list = []
bar = ProgBar(len(protein_fields), stream=sys.stdout, title='Collating per-protein lifestyle effects during regressions')
for protein in protein_fields:
    var_explained = read_csv('~/ch3_proteomics/3.1 covariates/regressing_out/effects/' + protein + '_effects_during_regressions.csv')
    vars_explained_list.append(var_explained)
    bar.update()
vars_explained_by_all_covars_at_each_step = concat(vars_explained_list)
vars_explained_by_all_covars_at_each_step.to_csv('~/ch3_proteomics/3.1 covariates/regressing_out/effects_during_regressions.csv', index=False)

# plot lifestyle effects over 'time'
vars_explained_by_all_covars_at_each_step = read_csv('~/ch3_proteomics/3.1 covariates/regressing_out/covariate_lifestyle_effects_during_regressions.csv')
vars_explained_by_all_covars_at_each_step['perc_var_explained'] = vars_explained_by_all_covars_at_each_step['var_explained'].apply(lambda x: x*100)

fig, ax = plt.subplots(figsize=[11.693, 8.268])
sns.lineplot(vars_explained_by_all_covars_at_each_step, x='stage', y='var_explained', hue='covar_name', estimator='mean', errorbar='sd',
             hue_order=['Well', 'Date at centre', 'Plate', 'Storage time', 'Date processed', 'Non-oily fish intake', 'Plate row', 'Plate column','Beef intake', 'Lamb intake', 'Pork intake', 'Water intake'])
ax.set_ylabel('R-squared value')
ax.set_title('R-squared value for models fitting lifestyles to data after each lifestyle is successively regressed out')
plt.legend(title='Covariate')
fig.tight_layout()
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/effects_during_regressions.png')
fig.savefig('~/ch3_proteomics/3.1 covariates/plots/effects_during_regressions.svg')
