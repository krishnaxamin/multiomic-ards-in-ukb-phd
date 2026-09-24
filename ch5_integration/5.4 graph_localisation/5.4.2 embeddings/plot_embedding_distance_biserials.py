""" Plot proportions of tissue-disease instances that have positive or negative rank biserials. """
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.use('TkAgg')

""" Read in data """
results_list = []
for dirr in [f"~/data/internal/knowledge_graph/gene_protein_distances/embeddings/{x.name}" for x in
             os.scandir('~/data/internal/knowledge_graph/gene_protein_distances/embeddings') if x.is_dir()]:
    results_list = results_list + [
        pd.read_csv(f"{dirr}/{x.name}").assign(proteomics=x.name.split('.')[0].split('_')[-1],
                                               gene_gene_links='_'.join(x.name.split('_')[2:-1])) for x in
        os.scandir(dirr) if x.is_file()]
results = pd.concat(results_list)

""" Collate data """
results_long = results[[x for x in results.columns if 'pval' not in x]].melt(
    id_vars=['disease', 'tissue', 'gwas_gene_disease_relationship', 'proteomics', 'gene_gene_links'],
    value_vars=[x for x in results.columns if 'inter_rank_biserial' in x], value_name='rank_biserial',
    var_name='comparison')
results_long['comparison'] = results_long['comparison'].replace(
    {'intra_gwas_vs_inter_rank_biserial': 'Intra-GWAS vs Inter',
     'intra_proteomics_vs_inter_rank_biserial': 'Intra-proteomics vs Inter'})
results_long_select = results_long[results_long['rank_biserial'].notna()]
results_long_select[
    'x_axis_category'] = results_long_select.gwas_gene_disease_relationship + '|' + results_long_select.proteomics + '|' + results_long_select.gene_gene_links
results_long_select['rank_biserial_binary'] = ['Positive' if x > 0 else 'Negative' for x in results_long_select.rank_biserial.to_list()]

total_cat_counts = results_long_select.value_counts(['x_axis_category', 'comparison']).reset_index()

summary_plot_df = pd.concat(
    [pd.DataFrame(sub_df.rank_biserial_binary.value_counts()).reset_index().assign(x_axis_category=cat[0], comparison=cat[1])
     for cat, sub_df in results_long_select.groupby(['x_axis_category', 'comparison'])])
summary_plot_df = summary_plot_df.merge(total_cat_counts.rename(columns={'count': 'total_count'}))
summary_plot_df['proportion'] = summary_plot_df['count'] * 100 / summary_plot_df['total_count']

# stats for text
summary_plot_df[(summary_plot_df['rank_biserial_binary'] == 'Positive') & (summary_plot_df['comparison'].str.contains('GWAS'))].proportion.describe()
summary_plot_df[(summary_plot_df['rank_biserial_binary'] == 'Positive') & (summary_plot_df['comparison'].str.contains('proteomics'))].proportion.describe()

""" Plot """
g = sns.FacetGrid(summary_plot_df, row='comparison', height=4, aspect=1.5, sharey=False)
g.map_dataframe(sns.barplot, x='x_axis_category', y='proportion', hue='rank_biserial_binary',
                palette=(plt.cm.coolwarm(1.0), plt.cm.coolwarm(0.0)),
                hue_order=['Positive', 'Negative'])
g.set_axis_labels('Options used in graph building (GWAS-gene links | proteomics algorithm | gene-gene links)',
                  'Proportion of \ntissue-ARD pairs')
g.set_titles(row_template='{row_name}')
g.tick_params(axis='x', rotation=90, labelsize=8)
for ax in g.axes_dict.values():
    ax.yaxis.grid(True, which='major', color='grey', alpha=0.6)
    ax.set_axisbelow(True)
g.add_legend()
g.tight_layout()
g.savefig(f"~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/rank_biserials.png")
g.savefig(f"~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/rank_biserials.svg")
plt.close()
