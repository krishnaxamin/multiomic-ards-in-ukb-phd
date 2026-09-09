""" Plot how many diseases genes are associated with via MAGMA, with spotlight on those with n_diseases >= 5 """
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.use('TkAgg')

adjustment = 'lifestyles'

magma_gene = pd.read_csv(f"~/data/internal/genomics/magma/{adjustment}/gene_level_results_assoc.csv")
magma_gene = magma_gene[magma_gene['bonf_sig'] == 1].copy()

disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(str) + ')'
# https://www.genenames.org/download/custom/
ncbi_to_symbol_mapping = pd.read_csv('https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit', sep='\t')[['Approved symbol', 'NCBI Gene ID']]
ncbi_to_symbol_mapping = ncbi_to_symbol_mapping[ncbi_to_symbol_mapping['NCBI Gene ID'].notna()].copy()

n_assoc_per_gene = pd.DataFrame(magma_gene.GENE.value_counts()).reset_index().set_axis(['GENE', 'gene_freq'], axis=1)

""" Plot n_assoc histogram """
plt.figure(figsize=(8, 6))
sns.histplot(n_assoc_per_gene['gene_freq'], log_scale=[None, 10], binwidth=1, discrete=True)
plt.xticks(list(range(1, n_assoc_per_gene.gene_freq.max() + 1)))
plt.grid(visible=True, axis='x', which='major', color='grey', alpha=0.2)
plt.xlabel('Number of diseases associated with each gene')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig(f"~/ch2_genomics/2.7 gene_level_analysis/plots/{adjustment}/n_diseases_per_gene.png")
plt.savefig(f"~/ch2_genomics/2.7 gene_level_analysis/plots/{adjustment}/n_diseases_per_gene.svg")
plt.close()

""" Plot the diseases with which genes that have multiple assoc diseases are associated """
magma_gene_focus = magma_gene[['GENE', 'disease']].merge(n_assoc_per_gene[n_assoc_per_gene.gene_freq >= 5]).merge(disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease').set_axis(['GENE', 'n_diseases', 'disease'], axis=1)
magma_gene_focus = magma_gene_focus.merge(ncbi_to_symbol_mapping.rename(columns={'Approved symbol': 'gene', 'NCBI Gene ID': 'GENE'})).drop('GENE', axis=1)
magma_gene_focus_heatmap = magma_gene_focus.assign(mark=1).pivot(index='gene', columns='disease', values='mark').fillna(0)

# Figure layout
fig = plt.figure(figsize=(8, 6))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)

# --- Scatter plot ---
ax = fig.add_subplot(gs[0])

sns.heatmap(
    data=magma_gene_focus_heatmap, linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
    cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
    linecolor='black'
)

ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Gene')
ax.set_xticklabels(ax.get_xticklabels(), size=8)
ax.set_yticklabels(ax.get_yticklabels(), size=8)

# --- Bar chart (counts per disease) ---
# log_vars = np.log1p(n_vars_per_disease_full.values)
magma_gene_focus_ax2 = magma_gene_focus[['n_diseases', 'gene']].drop_duplicates().sort_values(by='gene').reset_index(drop=True)
ax2 = fig.add_subplot(gs[1], sharey=ax)
ax2.barh(y=[x+0.5 for x in magma_gene_focus_ax2.index], width=magma_gene_focus_ax2['n_diseases'], color='grey', alpha=0.7)
ax2.set_xlabel('n(diseases)')
ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

# Add vertical grid lines at tick marks
ax2.set_xticks([2, 4, 6, 8])
ax2.grid(True, axis='x', color='lightgrey', lw=0.8, ls='--')

plt.savefig(f"~/ch2_genomics/2.7 gene_level_analysis/plots/{adjustment}/diseases_per_gene_with_many_diseases.png")
plt.savefig(f"~/ch2_genomics/2.7 gene_level_analysis/plots/{adjustment}/diseases_per_gene_with_many_diseases.svg")

plt.close()

""" Generate disease-gene table """
magma_gene_ncbi = magma_gene.merge(ncbi_to_symbol_mapping.rename(columns={'Approved symbol': 'gene', 'NCBI Gene ID': 'GENE'})).drop('GENE', axis=1)
magma_gene_ncbi = magma_gene_ncbi[['gene', 'disease']].copy()
magma_gene_ncbi = magma_gene_ncbi.merge(disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease')

magma_gene_ncbi_collated_list = []
for _, df in magma_gene_ncbi.groupby('icd10_three_letter'):
    magma_gene_ncbi_collated_list.append(df.assign(genes=', '.join(sorted(df.gene.to_list()))).drop(columns='gene').drop_duplicates())
magma_gene_ncbi_collated = pd.concat(magma_gene_ncbi_collated_list)
magma_gene_ncbi_collated = magma_gene_ncbi_collated[['icd10_three_letter', 'genes']].rename(columns={'icd10_three_letter': 'disease'})
magma_gene_ncbi_collated.to_csv(f"~/data/internal/genomics/magma/{adjustment}/disease_genes_mappings.csv", index=False)
