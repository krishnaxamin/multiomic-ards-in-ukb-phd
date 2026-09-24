""" Analyse co-occurrence (represented by alpha-MLE): distributions, Louvain clusters. """
import pandas as pd
import numpy as np
import scipy
import networkx as nx
import itertools
import upsetplot
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from statannotations.Annotator import Annotator
from utils.significance_labelling import significance_labelling

mpl.use('TkAgg')

""" Read in data """
genomics = pd.read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
proteomics = pd.read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
prior_proteomics = pd.read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences_alphamle.csv')

disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(str) + ')'
code_chapter_to_field_mapper = dict(zip(disease_info.code_chapter, disease_info.disease_field))

""" FDR data """
genomics = significance_labelling(genomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                  pair_columns=['disease1', 'disease2'])
genomics_sig = genomics[genomics['fdr_sig'] == 1].copy()
genomics_sig = genomics_sig.loc[~genomics_sig.index.duplicated(keep='first'), :]
proteomics = significance_labelling(proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                    pair_columns=['disease1', 'disease2'])
proteomics_sig = proteomics[proteomics['fdr_sig'] == 1].copy()
proteomics_sig = proteomics_sig.loc[~proteomics_sig.index.duplicated(keep='first'), :]
prior_proteomics = significance_labelling(prior_proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                          pair_columns=['disease1', 'disease2'])
prior_proteomics_sig = prior_proteomics[prior_proteomics['fdr_sig'] == 1].copy()
prior_proteomics_sig = prior_proteomics_sig.loc[~prior_proteomics_sig.index.duplicated(keep='first'), :]

""" Ranges """
genomics_sig.alpha_hat.describe()
proteomics_sig.alpha_hat.describe()
prior_proteomics_sig.alpha_hat.describe()

""" Plot alpha-MLE distributions """
all_alphas = pd.concat([genomics_sig[['alpha_hat']].assign(cat='Genomics'),
                        proteomics_sig[['alpha_hat']].assign(cat='Proteomics'),
                        prior_proteomics_sig[['alpha_hat']].assign(cat='Prior proteomics')])
pvalues = [
    scipy.stats.mannwhitneyu(genomics_sig.alpha_hat, proteomics_sig.alpha_hat).pvalue,
    scipy.stats.mannwhitneyu(genomics_sig.alpha_hat, prior_proteomics_sig.alpha_hat).pvalue,
    scipy.stats.mannwhitneyu(proteomics_sig.alpha_hat, prior_proteomics_sig.alpha_hat).pvalue
]
plt.figure(figsize=(6, 8))
ax = sns.violinplot(data=all_alphas, x='cat', y='alpha_hat', hue='cat', cut=0)
plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
plt.xlabel('Cohort')
plt.ylabel('Co-occurrence (alpha)')
# statistical annotations
pairs = [
    ('Genomics', 'Proteomics'),
    ('Genomics', 'Prior proteomics'),
    ('Proteomics', 'Prior proteomics')
]
annotator = Annotator(ax, pairs, data=all_alphas, x='cat', y='alpha_hat')
annotator.configure(test_short_name='MW', text_format='full', loc='outside')
# annotator._pvalue_format._pvalue_thresholds = [[1e-5, "1e-5"], [1e-4, "1e-4"],
#                                                [1e-3, "0.001"], [0.05 / 3, "0.017"]]
annotator.set_pvalues(pvalues)
annotator.annotate()
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.25)
# export
plt.tight_layout()
plt.savefig('~/ch4_other_data_layers/4.2 epidemiology/plots/alphamles_between_cohorts.png')
plt.savefig('~/ch4_other_data_layers/4.2 epidemiology/plots/alphamles_between_cohorts.svg')
plt.close()

""" Louvain """
# data_to_cluster = genomics_sig.copy()


def louvain(data_to_cluster: pd.DataFrame, omic: str):

    if omic == 'prior_proteomics':
        omic_dir = 'proteomics'
    else:
        omic_dir = omic

    data_for_graph = (data_to_cluster[['disease1', 'disease2', 'alpha_hat']].merge(
        disease_info[['disease_field', 'code_chapter']]
        .rename(columns={'disease_field': 'disease1'}))
                     .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                     .merge(disease_info[['disease_field', 'code_chapter']]
                            .rename(columns={'disease_field': 'disease2'}))
                     .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
    bma_graph = nx.Graph()

    # [(i, {'disease': <code_chapter>})]
    all_nodes = [(i, {'disease': x}) for i, x in enumerate(disease_info.code_chapter)]
    bma_graph.add_nodes_from(all_nodes)

    disease_node_mapper = {x[1]['disease']: x[0] for x in all_nodes}
    node_disease_mapper = {v: k for k, v in disease_node_mapper.items()}

    all_edges = []
    for i in range(len(data_for_graph)):
        disease1 = data_for_graph.at[i, 'disease1']
        disease2 = data_for_graph.at[i, 'disease2']
        alpha_hat = data_for_graph.at[i, 'alpha_hat']

        disease1_node = disease_node_mapper[disease1]
        disease2_node = disease_node_mapper[disease2]

        all_edges.append((disease1_node, disease2_node, {'alpha_hat': alpha_hat}))

    bma_graph.add_edges_from(all_edges)

    louvain_clusters = nx.algorithms.community.louvain_communities(bma_graph, seed=42, weight='alpha_hat')

    # take clusters and Louvain each one
    subclusters = []
    for cluster in louvain_clusters:
        cluster_graph = bma_graph.subgraph(cluster)
        louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='alpha_hat')
        subclusters = subclusters + louvain_subclusters

    # more Louvain - break down the clusters as far as they will go
    subclusters_next = []
    for cluster in subclusters:
        cluster_graph = bma_graph.subgraph(cluster)
        louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='alpha_hat')
        subclusters_next += louvain_subclusters

    while len(subclusters_next) != len(subclusters):
        subclusters = subclusters_next
        subclusters_next = []
        for cluster in subclusters:
            cluster_graph = bma_graph.subgraph(cluster)
            louvain_subclusters = nx.algorithms.community.louvain_communities(cluster_graph, seed=42, weight='alpha_hat')
            subclusters_next += louvain_subclusters

    # isolate non-trivial clusters
    louvain_groups = [[node_disease_mapper[x] for x in group] for group in [list(x) for x in subclusters_next if len(x) > 1]]

    # is multimorbidity higher within clusters than without - might this not suggest that multimorbidity is causal for related enrichments?
    multimorbidity = data_to_cluster.copy()

    multimorbidity = pd.merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
        columns={'disease_field': 'disease1'})).rename(columns={'disease1': 'disease1_field', 'code_chapter': 'disease1'})
    multimorbidity = pd.merge(multimorbidity, disease_info[['disease_field', 'code_chapter']].rename(
        columns={'disease_field': 'disease2'})).rename(columns={'disease2': 'disease2_field', 'code_chapter': 'disease2'})

    multimorbidity['disease_pair'] = multimorbidity.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                          axis=1)
    multimorbidity = multimorbidity.drop_duplicates('disease_pair')
    multimorbidity[['disease1', 'disease2']] = multimorbidity['disease_pair'].str.split('-', expand=True)


    # get multimorbidity scores for all disease pairs within each disease group
    louvain_groups_multimorbidity = {}
    cluster_combinations = []  # disease combinations which are present within clusters
    for group in louvain_groups:
        group_diseases_combinations = list(itertools.combinations(group, 2))
        group_alpha_hats = []
        for combo in group_diseases_combinations:
            multimorbidity_combo = multimorbidity.loc[
                (multimorbidity['disease1'] == sorted(combo)[0]) & (
                        multimorbidity['disease2'] == sorted(combo)[1]), 'alpha_hat']
            if len(multimorbidity_combo) == 0:
                continue
            alpha_hat = multimorbidity_combo.values[0]
            group_alpha_hats.append(alpha_hat)
        louvain_groups_multimorbidity['-'.join(group)] = {'alpha_hats': group_alpha_hats,
                                                          'mean_alpha_hat': np.mean(group_alpha_hats),
                                                          'median_alpha_hat': np.median(group_alpha_hats)}
        cluster_combinations = cluster_combinations + group_diseases_combinations

    # mean multimorbidity of whole dataset = 1.252
    mean_all_alpha_hat = multimorbidity['alpha_hat'].mean()

    # get alpha-hats for all pairs that are not in the same cluster (non-cluster pairs)
    cluster_combinations_idxs_alphamle = []
    for combo in cluster_combinations:
        multimorbidity_combo = multimorbidity.loc[
            (multimorbidity['disease1'] == sorted(combo)[0]) & (
                    multimorbidity['disease2'] == sorted(combo)[1]), 'alpha_hat']
        if len(multimorbidity_combo) == 0:
            continue
        cluster_combinations_idxs_alphamle.append(multimorbidity_combo.index.values[0])

    non_cluster_alpha_hats = multimorbidity.drop(cluster_combinations_idxs_alphamle).drop_duplicates('disease_pair')[
        'alpha_hat'].values

    louvain_groups_greater_alphamle = {}
    for group in louvain_groups:
        louvain_groups_greater_alphamle['-'.join(group)] = scipy.stats.mannwhitneyu(
            louvain_groups_multimorbidity['-'.join(group)]['alpha_hats'], non_cluster_alpha_hats,
            alternative='greater').pvalue

    """ Export Louvain groups """
    louvain_groups_as_fields = ['-'.join([code_chapter_to_field_mapper[x] for x in group]) for group in louvain_groups]
    louvain_group_df = pd.DataFrame.from_dict(louvain_groups_greater_alphamle, orient='index')
    louvain_group_df.columns = ['sig_greater_multimorbidity']
    louvain_group_df['bonf_sig_greater_multimorbidity'] = [1 if x < 0.05 / len(louvain_group_df) else 0
                                                           for x in louvain_group_df.sig_greater_multimorbidity.to_list()]
    louvain_group_df['disease_group'] = ['-'.join(sorted(group)) for group in louvain_groups]
    louvain_group_df[['disease_group', 'sig_greater_multimorbidity', 'bonf_sig_greater_multimorbidity']].to_csv(
        f"~/data/internal/{omic_dir}/{omic}_multimorbidity_disease_groups.csv", index=False)

    """ Return Louvain groups """
    louvain_group_df_export = louvain_group_df[louvain_group_df['bonf_sig_greater_multimorbidity'] == 1].reset_index(drop=True)
    louvain_group_df_export = louvain_group_df_export.assign(diseases=louvain_group_df_export['disease_group'].str.split('-'))
    louvain_group_df_export['disease_group'] = [int(x) + 1 for x in list(louvain_group_df_export.index)]

    return louvain_group_df_export[['disease_group', 'diseases']]


omic_dict = {'genomics': genomics_sig, 'proteomics': proteomics_sig, 'prior_proteomics': prior_proteomics_sig}
omic_results_dict = {}
for omic, omic_df in omic_dict.items():
    omic_results_dict[omic] = louvain(data_to_cluster=omic_df, omic=omic)

""" UpSet """


def scaled_upset_plot(upsetplot_data):
    # given upsetdata (generated from from_contents or similar), generate parameters that scale the upsetplot
    #  proportionally to visually fit everything properly

    # upsetplot_data will have index of form (True, True, False, ...)
    #  this gets the length of that boolean vector = number of sets there are
    n_sets = len(upsetplot_data.index.names)
    n_intersections = len(set(upsetplot_data.index))

    # figure size
    fig_width = np.clip(n_intersections * 0.6, 8, 40)
    fig_height = np.clip(n_sets * 0.5, 5, 20)

    # text size: shrinks as sets increase, floored at 5pt
    label_size = np.clip(12 - (n_sets - 5) * 0.3, 5, 12)

    # dot and line size: shrink as the matrix gets denser
    dot_size = np.clip(200 - (n_sets - 5) * 6, 20, 200)
    linewidth = np.clip(1.2 - (n_sets - 5) * 0.03, 0.4, 1.2)

    omics_colours = {'genomics': '#1f77b4', 'proteomics': '#ff7f0e', 'prior_proteomics': '#2ca02c'}

    fig = plt.figure(figsize=(fig_width, fig_height))
    gene_set_dict_working_upsetplot = upsetplot.UpSet(upsetplot_data, element_size=None, show_counts=True,
                                                      sort_categories_by='input')
    # gene_set_dict_working_upsetplot.style_subsets(min_degree=n_sets,
    #                                               facecolor='#1f77b4')  # colour by 1st tab10 colour
    # colour the totals bar charts
    for omic, colour in omics_colours.items():
        gene_set_dict_working_upsetplot.style_categories(categories=[x for x in list(upsetplot_data.index.names)
                                                                     if x.split('-')[0] == omic],
                                                         bar_facecolor=colour)

    upsetplot_axes_dict = gene_set_dict_working_upsetplot.plot(fig=fig)

    ax_matrix = upsetplot_axes_dict['matrix']
    for collection in ax_matrix.collections:
        if isinstance(collection, mpl.collections.PathCollection):
            collection.set_sizes([dot_size])
        elif isinstance(collection, mpl.collections.LineCollection):
            collection.set_linewidths(linewidth)

        # colour dots by category (genomics, proteomics, prior_proteomics)
        offsets = collection.get_offsets()
        facecolors = collection.get_facecolor()

        # Skip empty collections (e.g. connector lines)
        if len(facecolors) == 0 or len(offsets) == 0:
            continue

        new_colors = []
        for i, (x, y) in enumerate(offsets):
            # Check if dot is already filled (not grey/white background)
            existing = facecolors[i] if len(facecolors) > 1 else facecolors[0]
            is_filled = existing[3] > 0.5  # alpha check
            if is_filled:
                cat_index = int(round(y))
                cat_name = upsetplot_data.index.names[cat_index]
                new_colors.append(omics_colours.get(cat_name.split('-')[0], 'black'))
            else:
                new_colors.append(existing)  # keep original unfilled style
        collection.set_facecolor(new_colors)

    # for collection in ax_matrix.collections:
    #     collection.set_sizes([dot_size])      # scatter circles
    # for line in ax_matrix.lines:
    #     line.set_linewidth(linewidth)         # connecting lines

    ax_matrix.tick_params(axis='y', labelsize=label_size)
    upsetplot_axes_dict['intersections'].tick_params(labelsize=label_size * 0.9)
    if upsetplot_axes_dict.get('totals') is not None:
        upsetplot_axes_dict['totals'].tick_params(labelsize=label_size * 0.9)

    legend_handles = [
        mpl.patches.Patch(color=color, label=cat)
        for cat, color in omics_colours.items()
    ]
    ax_matrix.legend(handles=legend_handles, bbox_to_anchor=(1.04, 0.5), loc='center left', title='Cohort')

    fig.tight_layout()
    return fig, upsetplot_axes_dict


groupings_upset = upsetplot.from_contents(dict(zip([f"genomics-group{x}" for x in list(omic_results_dict['genomics'].disease_group.to_list())],
                                                   omic_results_dict['genomics'].diseases.to_list())) |
                                          dict(zip([f"proteomics-group{x}" for x in list(omic_results_dict['proteomics'].disease_group.to_list())],
                                                   omic_results_dict['proteomics'].diseases.to_list())) |
                                          dict(zip([f"prior_proteomics-group{x}" for x in list(omic_results_dict['prior_proteomics'].disease_group.to_list())],
                                                   omic_results_dict['prior_proteomics'].diseases.to_list())))

groupings_upset_upsetplot, axes = scaled_upset_plot(groupings_upset)
groupings_upset_upsetplot.savefig(f"~/ch4_other_data_layers/4.2 epidemiology/plots/louvain_groupings_upset.svg", bbox_inches='tight')
groupings_upset_upsetplot.savefig(f"~/ch4_other_data_layers/4.2 epidemiology/plots/louvain_groupings_upset.png", bbox_inches='tight')
plt.close()
