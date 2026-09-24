"""
Graph to show the number of FDR < 0.05 metabolites per disease, and a heatmap of their effect sizes.
"""
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('TkAgg')

""" Read in data """
metabolomics = pd.read_csv(f"~/data/internal/metabolomics/pan-ukbb-eur_assoc_metabolomics_firth.csv")
sig = metabolomics[metabolomics.fdr_sig == 1].copy()

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

sig = sig.merge(disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease')
sig = sig[['code_chapter', 'metabolite', 'beta']].copy()

"""" Get metabolite counts """
n_metabolites_per_disease = pd.DataFrame(sig.value_counts('code_chapter')).reset_index().sort_values(by='code_chapter').reset_index(drop=True)

""" Plot n(metabolites) per disease """
# Figure layout
fig = plt.figure(figsize=(12, 14))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], wspace=0.05)

# --- Heatmap ---
ax = fig.add_subplot(gs[1])

sig_heatmap = sig.pivot_table(values='beta', index='metabolite', columns='code_chapter').fillna(0)

sns.heatmap(sig_heatmap, cmap='coolwarm', center=0, yticklabels=False, xticklabels=1, ax=ax,
            cbar_kws={'label': 'Beta'})
ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Metabolite')

# --- Bar chart (metabolite counts per disease) ---
n_metabolites_per_disease['log_count'] = np.log1p(n_metabolites_per_disease['count'])
# n_metabolites_per_disease['code_chapter_pos'] = np.arange(len(n_metabolites_per_disease)) + 0.5
ax2 = fig.add_subplot(gs[0], sharex=ax)
ax2.bar(np.arange(len(n_metabolites_per_disease)) + 0.5, n_metabolites_per_disease['log_count'], color='grey', alpha=0.7)
ax2.set_ylabel('ln(n_metabolites + 1)')
ax2.tick_params(axis='x', rotation=90)

# Add raw counts labels
for x, count in zip(range(len(n_metabolites_per_disease)), n_metabolites_per_disease['count'].values):
    ax2.text(x + 0.5, n_metabolites_per_disease['log_count'][x] + 0.05, str(count), ha='center', fontsize=5)

plt.show()

plt.savefig(f"~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/plots/assoc_metabolites_all_disease_firth.png")
plt.savefig(f"~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/plots/assoc_metabolites_all_disease_firth.svg")
plt.close()

""" Get n(diseases) per metabolite """
n_assoc_per_metabolite = pd.DataFrame(sig.metabolite.value_counts()).reset_index().set_axis(['metabolite', 'n_diseases'], axis=1)

""" Plot n_assoc histogram """
plt.figure(figsize=(8, 6))
sns.histplot(n_assoc_per_metabolite.n_diseases, log_scale=[None, 10], binwidth=1, discrete=True)
plt.xticks(list(range(1, n_assoc_per_metabolite.n_diseases.max() + 1)))
plt.grid(visible=True, axis='x', which='major', color='grey', alpha=0.2)
plt.xlabel('Number of diseases associated with each metabolite')
plt.ylabel('Frequency')
plt.gca().set_xticks([1] + list(range(5, math.ceil(n_assoc_per_metabolite.n_diseases.max()/5 + 1) * 5, 5)))
plt.tight_layout()
plt.savefig(f"~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/plots/n_diseases_per_metabolite_firth.png")
plt.savefig(f"~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/plots/n_diseases_per_metabolite_firth.svg")
plt.close()
