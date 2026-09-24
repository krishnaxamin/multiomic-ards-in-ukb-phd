""" Plot the number of clusters for each tissue/gene-gene link combination. """
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.use('TkAgg')

""" Read in data """
results_list = []
for dirr in [f"~/data/internal/knowledge_graph/tissue_graphs/{x.name}" for x in os.scandir('~/data/internal/knowledge_graph/tissue_graphs') if x.is_dir()]:
    results_list = results_list + [pd.read_csv(f"{dirr}/{x.name}").assign(gene_gene_links='_'.join(x.name.split('_')[1:-1]),
                                                                          tissue=os.path.basename(dirr))
                                   for x in os.scandir(dirr) if x.is_file() and '_firth.csv' in x.name]
results = pd.concat(results_list)

""" Count non-trivial clusters per gene-gene link and tissue """
results_non_trivial = results[results['cluster_genes'].str.contains('|', regex=False)].copy()

n_clusters = results_non_trivial.value_counts(['tissue', 'gene_gene_links']).reset_index().set_axis(['tissue', 'gene_gene_link', 'counts'], axis=1)
n_clusters.counts.describe()

""" Plot counts """
n_clusters['tissue'] = pd.Categorical(n_clusters['tissue'], sorted(list(set(n_clusters.tissue))))
g = sns.FacetGrid(n_clusters, row='gene_gene_link', aspect=2.5, sharey=False)
g.map_dataframe(sns.barplot, x='tissue', y='counts')
g.set_axis_labels('Tissue', 'Number of clusters')
g.set_titles(row_template='{row_name}')
g.tick_params(axis='x', rotation=90)
for ax in g.axes_dict.values():
    ax.yaxis.grid(True, which='major', color='grey', alpha=0.6)
    ax.set_axisbelow(True)
g.tight_layout()
g.savefig(f"~/ch5_integration/5.4 graph_localisation/5.4.1 clusters/plots/n_clusters_per_tissue.png")
g.savefig(f"~/ch5_integration/5.4 graph_localisation/5.4.1 clusters/plots/n_clusters_per_tissue.svg")
plt.close()
