"""
Plot HDL and LDSC correlation results for lifestyles-adjusted GWAS on a shared heatmap
"""
from utils.significance_labelling import significance_labelling

import pandas as pd
import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('TkAgg')

plot_path = '~/ch2_genomics/2.10 whole_genome_genetic_correlations/hdl_ldsc'

""" Read in results """

# HDL
hdl = pd.read_csv('~/data/internal/genomics/whole_genome_genetic_correlations/hdl_genetic_correlations.csv')

# LDSC
ldsc = pd.read_csv('~/data/internal/genomics/whole_genome_genetic_correlations/ldsc_genetic_correlations.csv')
ldsc.loc[ldsc['genetic_correlation'] > 1, 'genetic_correlation'] = 1
ldsc.loc[ldsc['genetic_correlation'] < -1, 'genetic_correlation'] = -1

# disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Significance label """
hdl = significance_labelling(hdl, pvalue_col='pval')
ldsc = significance_labelling(ldsc, pvalue_col='pval')

""" Prep for plotting """
# top LH half = HDL (disease1 < disease2)
# bottom RH half = LDSC (disease1 > disease2)

data_to_plot = pd.concat([hdl, ldsc.rename(columns={'disease1': 'disease2', 'disease2': 'disease1'})])[
    ['disease1', 'disease2', 'genetic_correlation', 'fdr_sig']]
data_to_plot = (data_to_plot
                .merge(disease_info[['disease_field', 'code_chapter']]
                       .rename(columns={'disease_field': 'disease1'}))
                .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                .merge(disease_info[['disease_field', 'code_chapter']]
                       .rename(columns={'disease_field': 'disease2'}))
                .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                          categories=sorted(data_to_plot['disease1'].unique()))
data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                          categories=sorted(data_to_plot['disease2'].unique()))
data_to_plot['fdr_sig'] = data_to_plot['fdr_sig'].astype(str)
data_to_plot.loc[data_to_plot['fdr_sig'] == '1', 'fdr_sig'] = 'Yes'
data_to_plot.loc[data_to_plot['fdr_sig'] == '0', 'fdr_sig'] = 'No'

""" Plot """
# fill = genetic_correlation, shape = fdr_sig

# map data to symmetric categorical index space
all_labels = sorted(set(data_to_plot['disease1']) | set(data_to_plot['disease2']))
# make mapping have 1-tick gaps between chapters, so that there are whitespace gaps/boundaries between chapters
label_to_num = {}
current_pos = 0
for i in range(len(all_labels)):
    if i == 0:
        label_to_num[all_labels[i]] = current_pos
        current_pos += 1
        continue
    current_chapter = int(re.search(r"\((\d+)\)", all_labels[i]).group(1))
    past_chapter = int(re.search(r"\((\d+)\)", all_labels[i - 1]).group(1))
    if past_chapter < current_chapter:
        current_pos += 1  # insert a 1-tick wide gap into the axis at the boundary between chapters
    label_to_num[all_labels[i]] = current_pos
    current_pos += 1
# label_to_num = {lab: i for i, lab in enumerate(all_labels)}
data_to_plot['disease1_as_num'] = data_to_plot['disease1'].map(label_to_num)
data_to_plot['disease2_as_num'] = data_to_plot['disease2'].map(label_to_num)

# plot
fig, ax = plt.subplots(figsize=(12, 10))
sns.scatterplot(data=data_to_plot.rename(
    columns={'fdr_sig': 'FDR < 0.05', 'genetic_correlation': 'Genetic correlation'}
),
    x='disease1_as_num', y='disease2_as_num', hue='Genetic correlation', style='FDR < 0.05',
    style_order=['Yes', 'No'], zorder=3, palette=plt.cm.coolwarm, hue_norm=(-1, 1),
    edgecolor=None, markers={'Yes': 'o', 'No': '^'}, ax=ax)
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.grid(True, zorder=0)

if ax.yaxis_inverted():
    ax.invert_yaxis()

# set ticks labels to the string versions
tick_positions = [label_to_num[label] for label in all_labels]
ax.set_xticks(tick_positions)
ax.set_xticklabels(all_labels, rotation=90, size=9)
ax.set_yticks(tick_positions)
ax.set_yticklabels(all_labels, size=9)

# add diagonal line going top-left to bottom-right, ensuring whitespaces are dealt with properly
coords = np.array(tick_positions)
ax.plot(coords, coords,
        color='black', linestyle='--', linewidth=0.5)

# bold legend section titles
legend = ax.get_legend()
for text, handle in zip(legend.texts, legend.legend_handles):
    if handle._label in ['FDR < 0.05', 'Genetic correlation']:  # section header
        text.set_fontweight('bold')

plt.tight_layout()
plt.savefig(f"{plot_path}.png")
plt.savefig(f"{plot_path}.svg")
plt.close()
