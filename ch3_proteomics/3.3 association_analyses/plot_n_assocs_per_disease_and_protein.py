"""
Graph to show the number of FDR < 0.05 proteins per disease, and a heatmap of their effect sizes.
"""
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('TkAgg')

algorithm = 'cox'

""" Read in data """
proteomics = pd.read_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{algorithm}.csv")
sig = proteomics[proteomics.fdr_sig == 1].copy()

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

sig = sig.merge(disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease')
sig = sig[['code_chapter', 'protein', 'beta']].copy()

"""" Get protein counts """
n_proteins_per_disease = pd.DataFrame(sig.value_counts('code_chapter')).reset_index().sort_values(by='code_chapter').reset_index(drop=True)

""" Plot n(proteins) per disease """
# Figure layout
fig = plt.figure(figsize=(12, 14))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], wspace=0.05)

# --- Heatmap ---
ax = fig.add_subplot(gs[1])

sig_heatmap = sig.pivot_table(values='beta', index='protein', columns='code_chapter').fillna(0)

sns.heatmap(sig_heatmap, cmap='coolwarm', center=0, yticklabels=False, xticklabels=1, ax=ax,
            cbar_kws={'label': 'Beta'})
ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Protein')

# --- Bar chart (protein counts per disease) ---
n_proteins_per_disease['log_count'] = np.log1p(n_proteins_per_disease['count'])
# n_proteins_per_disease['code_chapter_pos'] = np.arange(len(n_proteins_per_disease)) + 0.5
ax2 = fig.add_subplot(gs[0], sharex=ax)
ax2.bar(np.arange(len(n_proteins_per_disease)) + 0.5, n_proteins_per_disease['log_count'], color='grey', alpha=0.7)
ax2.set_ylabel('ln(n_proteins + 1)')
ax2.tick_params(axis='x', rotation=90)

# Add raw counts labels
for x, count in zip(range(len(n_proteins_per_disease)), n_proteins_per_disease['count'].values):
    if algorithm == 'cox':
        ax2.text(x + 0.5, n_proteins_per_disease['log_count'][x] + 0.05, str(count), ha='center', fontsize=5)
    else:
        ax2.text(x + 0.5, n_proteins_per_disease['log_count'][x] + 0.05, str(count), ha='center', fontsize=8)

plt.show()

plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/assoc_proteins_all_disease_{algorithm}.png")
plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/assoc_proteins_all_disease_{algorithm}.svg")
plt.close()

""" Get n(diseases) per protein """
n_assoc_per_protein = pd.DataFrame(sig.protein.value_counts()).reset_index().set_axis(['protein', 'n_diseases'], axis=1)

""" Plot n_assoc histogram """
plt.figure(figsize=(8, 6))
sns.histplot(n_assoc_per_protein.n_diseases, log_scale=[None, 10], binwidth=1, discrete=True)
plt.xticks(list(range(1, n_assoc_per_protein.n_diseases.max() + 1)))
plt.grid(visible=True, axis='x', which='major', color='grey', alpha=0.2)
plt.xlabel('Number of diseases associated with each protein')
plt.ylabel('Frequency')
if algorithm == 'cox':
    plt.gca().set_xticks([1] + list(range(5, math.ceil(n_assoc_per_protein.n_diseases.max()/5 + 1) * 5, 5)))
plt.tight_layout()
plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/n_diseases_per_protein_{algorithm}.png")
plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/n_diseases_per_protein_{algorithm}.svg")
plt.close()

""" Plot the diseases with which proteins that have 5+ assoc diseases are associated """
if algorithm == 'firth':
    n_disease_limit = 5
    x_tick_fontsize = 10
    y_tick_fontsize = 10
elif algorithm == 'cox':
    n_disease_limit = 45
    x_tick_fontsize = 8
    y_tick_fontsize = 9
else:
    raise ValueError('Wrong algorithm')
disease_protein_focus = sig[['protein', 'code_chapter']].merge(n_assoc_per_protein[n_assoc_per_protein.n_diseases >= n_disease_limit])
disease_protein_focus_heatmap = disease_protein_focus.assign(mark=1).pivot(index='protein', columns='code_chapter', values='mark').fillna(0)

# Figure layout
fig = plt.figure(figsize=(8, 6))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)

# --- Scatter plot ---
ax = fig.add_subplot(gs[0])

sns.heatmap(
    data=disease_protein_focus_heatmap, linewidths=0.15, ax=ax, xticklabels=1, yticklabels=1,
    cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
    linecolor='black'
)

ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Protein')
# ax.set_xticklabels(ax.get_xticklabels(), size=8 - math.floor(disease_protein_focus.code_chapter.nunique() / 50))
ax.set_xticklabels(ax.get_xticklabels(), size=x_tick_fontsize)
ax.set_yticklabels(ax.get_yticklabels(), size=y_tick_fontsize)

# --- Bar chart (counts per disease) ---
disease_protein_focus_ax2 = disease_protein_focus[['n_diseases', 'protein']].drop_duplicates().sort_values(by='protein').reset_index(drop=True)
ax2 = fig.add_subplot(gs[1], sharey=ax)
ax2.barh(y=[x+0.5 for x in disease_protein_focus_ax2.index], width=disease_protein_focus_ax2['n_diseases'], color='grey', alpha=0.7)
ax2.set_xlabel('n(diseases)')
ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

# Add vertical grid lines at tick marks
ax2.grid(True, which='major', axis='x', color='lightgrey', lw=0.8, ls='--')
ax2.set_axisbelow(True)

# Add raw counts labels
for x, count in zip(range(len(disease_protein_focus_ax2)), disease_protein_focus_ax2['n_diseases'].values):
    ax2.text(count + 0.05, x + 0.5, str(count), va='center', fontsize=8)
ax2.set_xlim(ax2.get_xlim()[0], ax2.get_xlim()[1] + (1*math.floor(ax2.get_xlim()[1]/15)))

plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/diseases_per_protein_with_many_diseases_{algorithm}.png")
plt.savefig(f"~/ch3_proteomics/3.3 association_analyses/plots/diseases_per_protein_with_many_diseases_{algorithm}.svg")

plt.close()
