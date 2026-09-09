""" Construct LCV results as a DAG and calculate centrality metrics for each ARD/node from it. """
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

mpl.use('TkAgg')

lcv = pd.read_csv('~/data/internal/genomics/causality/lcv_results_processed.csv')
lcv = lcv[(abs(lcv.gcp_est) > 0.6) & (lcv.gcp0_fdr < 0.05)].copy()

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

lcv = (lcv
       .merge(disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease1'}))
       .drop(columns='disease1').rename(columns={'icd10_three_letter': 'disease1'})
       .merge(disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease2'}))
       .drop(columns='disease2').rename(columns={'icd10_three_letter': 'disease2'}))[['disease1', 'disease2', 'gcp_est']].copy()

""" Create graph """
# initialise graph
G = nx.DiGraph()

# add nodes and edges
for row in lcv.to_dict(orient='records'):
    if row['gcp_est'] < 0:
        src = row['disease2']
        tgt = row['disease1']
    else:
        src = row['disease1']
        tgt = row['disease2']

    G.add_edge(src, tgt)

""" Compute centralities """
in_degree = dict(G.in_degree())
out_degree = dict(G.out_degree())

# Betweenness on a disconnected graph is fine; nx handles it per-component.
betweenness = nx.betweenness_centrality(G, normalized=True)

# PageRank – use a small alpha for sparse graphs; 0.85 is the classic default.
try:
    pagerank = nx.pagerank(G, alpha=0.85, max_iter=1000)
except nx.PowerIterationFailedConvergence:
    print('Warning: PageRank did not converge – results may be approximate.')
    pagerank = nx.pagerank(G, alpha=0.85, max_iter=1000, tol=1e-3)

rows = []
for node in sorted(G.nodes(), key=str):
    rows.append({
        'node': node,
        'in_degree': in_degree[node],
        'out_degree': out_degree[node],
        'betweenness': round(betweenness[node], 6),
        'pagerank': round(pagerank[node], 6),
    })
centrality_df = pd.DataFrame(rows)
centrality_df.to_csv('~/data/internal/genomics/causality/dag_centrality.csv', index=False)

""" Plot centrality metrics """
centrality_to_plot = centrality_df.merge(disease_info[['icd10_three_letter', 'code_chapter']].rename(columns={'icd10_three_letter': 'node'})).drop(columns='node')
centrality_to_plot['code_chapter'] = pd.Categorical(centrality_to_plot.code_chapter, sorted(list(centrality_to_plot.code_chapter)))

metrics = ['in_degree', 'out_degree', 'betweenness', 'pagerank']
metrics_names = ['In-degree', 'Out-degree', 'Betweenness', 'PageRank']

for col in metrics:
    highlighting = centrality_to_plot[['code_chapter', col]].sort_values(by=col, ascending=False).iloc[:3, :]
    highlighting[f"{col}_highlight"] = 1
    centrality_to_plot = centrality_to_plot.merge(highlighting, how='left')
    centrality_to_plot[f"{col}_highlight"] = centrality_to_plot[f"{col}_highlight"].fillna(0)

fig = plt.figure(figsize=(6, 6))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1], hspace=0.05)

ax = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1], sharex=ax)
ax2 = fig.add_subplot(gs[2], sharex=ax)
ax3 = fig.add_subplot(gs[3], sharex=ax)
axes = [ax, ax1, ax2, ax3]

for i in range(len(metrics)):
    sns.barplot(data=centrality_to_plot, x='code_chapter', y=metrics[i], hue=f"{metrics[i]}_highlight",
                hue_order=[0, 1], palette=['grey', plt.cm.tab10(i)], ax=axes[i])
    if i == len(metrics) - 1:
        axes[i].set_xlabel('Disease ICD-10 code (ICD-10 chapter)', fontsize=10)
        axes[i].tick_params(axis='x', rotation=90, labelsize=12)
    else:
        axes[i].set_xlabel('')
        axes[i].tick_params(axis='x', labelbottom=False)
    axes[i].tick_params(axis='y', labelsize=10)
    axes[i].get_legend().remove()
    axes[i].set_ylabel(metrics_names[i], fontsize=12)

patches = [mpl.patches.Patch(color=plt.cm.tab10(i), label=metrics_names[i]) for i in range(len(metrics))]

fig.savefig('~/ch2_genomics/2.14 causality/plots/dag_centrality.png')
fig.savefig('~/ch2_genomics/2.14 causality/plots/dag_centrality.svg')
plt.close()
