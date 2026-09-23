"""
Graph to show the number of Bonferroni-adjusted significant INFO > 0.8 variants and genes (via MAGMA) per disease,
as well as the distribution of those associated variants across the genome.
Requires completion of MAGMA analysis and results processing.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import numpy as np

mpl.use('TkAgg')

adjustment = 'minimal'

""" Read in variant data """
variants = pd.read_csv(f"~/data/internal/genomics/pan-ukbb-eur_assoc_saige_{adjustment}.csv", dtype={'CHR': str, 'POS': int})
variants = variants[(variants['assoc_label'] == 'assoc') & (variants['info_label'] == 'high_info')][['CHR', 'POS', 'MarkerID', 'disease']]

""" Read in chromosome hg19 coordinate data """
chr_coords = pd.read_csv('~/data/external/genomics/chr_coordinates_hg19.csv')
chr_coords['midpoint'] = (chr_coords['from'] + chr_coords['to']) / 2
chr_coords['chr_num'] = chr_coords['chr'].apply(lambda x: x.replace('chr', ''))
chr_offsets = dict(zip(list(chr_coords['chr_num']), list(chr_coords['from'])))

""" Prep variant data for plotting """
variants['pos_cum'] = variants.apply(lambda r: r['POS'] + chr_offsets[r['CHR']], axis=1)

n_vars_per_disease = variants.groupby('disease').size()

""" Read in MAGMA gene data """
magma_results = pd.read_csv('~/data/internal/genomics/magma/' + adjustment + '/gene_level_results_assoc.csv')
num_assocs_per_disease = pd.concat(
    [magma_results.groupby(['disease']).nom_sig.sum().to_frame(name='nom_sig'),
     magma_results.groupby(['disease']).fdr_sig.sum().to_frame(name='fdr_sig'),
     magma_results.groupby(['disease']).bonf_sig.sum().to_frame(name='bonf_sig'),
     magma_results.groupby(['disease']).fdr_sig_permissive.sum().to_frame(name='fdr_sig_permissive')],
    axis=1)
num_assocs_per_disease['disease'] = num_assocs_per_disease.index

disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
magma_results = magma_results.merge(disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease').rename(columns={'icd10_three_letter': 'disease'})
n_genes_per_disease = pd.Series(magma_results['bonf_sig'].values, index=magma_results['disease'])
n_genes_per_disease = n_genes_per_disease[n_genes_per_disease != 0]

""" Align which diseases for which to plot data """
# all diseases with non-zero data in one of the datasets
all_diseases_with_data = pd.Index(
    np.unique(
        np.concatenate([
            variants['disease'].unique(),         # variants scatter
            n_vars_per_disease.index.values,    # variant counts
            n_genes_per_disease.index.values        # gene counts
        ])
    )
)

""" Set order in which diseases are plotted """
disease_order = sorted(list(all_diseases_with_data))

# reindex all series datasets
n_vars_per_disease_full = n_vars_per_disease.reindex(all_diseases_with_data, fill_value=0)
n_genes_per_disease_full = n_genes_per_disease.reindex(all_diseases_with_data, fill_value=0)

""" Plot """
# Figure layout
fig = plt.figure(figsize=(14, 12))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 3, width_ratios=[4, 1, 1], wspace=0.05)

# --- Scatter plot ---
ax = fig.add_subplot(gs[0])

sns.stripplot(
    data=variants, x='pos_cum', y='disease',
    order=disease_order,
    jitter=0.3, size=3, alpha=0.6, color='black', ax=ax
)

# Chromosome alternating background, boundary lines
for i, row in chr_coords.iterrows():
    # Alternate background color for chromosomes
    ax.axvspan(row['from'], row['to'], color='lightgrey' if i % 2 == 0 else 'white', alpha=0.3, zorder=-1)

    # Vertical line at chromosome boundary (right edge)
    ax.axvline(x=row['to'], color='grey', lw=0.5, ls='--', alpha=0.7)

# Remove normal x ticks and x-axis label
ax.set_xticks([])
ax.set_xlabel('')

ax.set_ylabel('Disease')

# --- Bar chart (counts per disease) ---
log_vars = np.log1p(n_vars_per_disease_full.values)
ax2 = fig.add_subplot(gs[1], sharey=ax)
ax2.barh(n_vars_per_disease_full.index, log_vars, color='grey', alpha=0.7)
ax2.set_xlabel('ln(n_vars + 1)')
ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

# Add vertical grid lines at tick marks
xticks = ax2.get_xticks()
for xt in xticks:
    ax2.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)
    
# Add raw counts labels
for y, count in zip(range(len(n_vars_per_disease_full)), n_vars_per_disease_full.values):
    ax2.text(log_vars[y] + 0.05, y, str(count), va='center', fontsize=8)

# --- Bar chart (gene counts per disease) ---
log_genes = np.log1p(n_genes_per_disease_full.values)
ax3 = fig.add_subplot(gs[2], sharey=ax)
ax3.barh(n_genes_per_disease_full.index, log_genes, color='darkgrey', alpha=0.7)
ax3.set_xlabel('ln(n_genes + 1)')
ax3.tick_params(axis='y', left=False, labelleft=False)

# Vertical grid lines
for xt in ax3.get_xticks():
    ax3.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

# Add raw counts labels
for y, count in zip(range(len(n_genes_per_disease_full)), n_genes_per_disease_full.values):
    ax3.text(log_genes[y] + 0.05, y, str(count), va='center', fontsize=8)

# --- Add chromosome labels on the variant distribution section ------
for _, row in chr_coords.iterrows():
    ax.text(
        row['midpoint'],               # centered under chromosome
        ax.get_ylim()[0] + 0.2,        # just below the bottom of the plot
        str(row['chr_num']),
        ha='center', va='top',
        fontsize=8, color='black',
        clip_on=False
    )

# Add x-axis label
ax.text(
    0.5, -0.025, 'Chromosome',
    ha='center', va='top', transform=ax.transAxes, fontsize=10
)

plt.show()

plt.savefig(f"~/ch2_genomics/2.4 process_gwas_results/plots/assoc_vars_all_disease_{adjustment}.png")
plt.savefig(f"~/ch2_genomics/2.4 process_gwas_results/plots/assoc_vars_all_disease_{adjustment}.svg")
plt.close()
