from utils.enrichment_analyses import prep_annotation_files
from statannotations.Annotator import Annotator
from concepts import Context
from pyprind import ProgBar

import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
import pandas as pd
import itertools
import math


mpl.use('TkAgg')

""" Read invariant data """

# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
disease_info_slim = disease_info[['disease_field', 'icd10_three_letter', 'icd10_chapter']].copy()
code_chapter_to_field_mapper = dict(zip(disease_info['code_chapter'], disease_info['disease_field']))
field_to_code_chapter_mapper = dict(zip(disease_info['disease_field'], disease_info['code_chapter']))

# ontology data
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/NCBI2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
annotations_ic = -np.log(annotations.annotation_id.value_counts() / len(annotations.id.unique()))

""" Read in data """
bma_mica_contributions = pd.read_csv(f"~/data/internal/genomics/pathway_based_similarities/disease_resnik_micas_reactome.csv")

""" Functions """


def count_high_flying_micas(selected_micas_df: pd.DataFrame) -> pd.DataFrame:
    high_flying_micas_slim = selected_micas_df[['mica', 'disease1', 'disease2']].drop_duplicates()
    high_flying_micas_diseases = [
        {'mica': mica,
         'diseases': '|'.join([field_to_code_chapter_mapper[x]
                               for x in
                               sorted(list(set(mica_diseases_df['disease1']) | set(mica_diseases_df['disease2'])))])}
        for mica, mica_diseases_df in high_flying_micas_slim.groupby('mica')]
    high_flying_micas_diseases_df = pd.DataFrame(high_flying_micas_diseases)
    high_flying_micas_diseases_df['n_diseases'] = high_flying_micas_diseases_df['diseases'].apply(
        lambda x: len(x.split('|')))
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        pd.DataFrame(high_flying_micas_slim['mica'].value_counts().rename('pair_count')).reset_index())

    # get reactome pathway names
    reactome_pathways = pd.read_csv('graph_construction/database_data_files/reactome/ReactomePathways.txt', sep='\t',
                                    header=None, names=['id', 'name', 'species'])
    reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

    return high_flying_micas_diseases_df.sort_values(by='pair_count', ascending=False).reset_index(drop=True)


def plot_high_flying_micas(counted_micas_df: pd.DataFrame, n_top_pathways: int, plot_name: str, plot_dir: str):
    """
    Given a set of high-flying driver MICAs, plot their links with diseases and disease pairs via a
    (1) heatmap (MICA-disease)
    (2) bar chart counting n(diseases)
    (3) bar chart counting n(disease pairs)
    The MICAs are ordered by n(disease pairs) because we are interested in MICAs that drive similarity between pairs.
    :param plot_dir:
    :param counted_micas_df:
    :param n_top_pathways:
    :param plot_name:
    :return:
    """
    # heatmap. mica_name on y-axis (cut-off after n characters), ARDs x-axis.
    # Separate bar charts of n_diseases and pair_count, as with the combined manhattan-ish plot
    counted_micas_df = counted_micas_df.reset_index(drop=True)

    # split mica_name over 2 lines if n_characters > 50
    def mica_name_split(name):
        if len(name) <= 50:
            return name

        lines = []
        current_line = ""
        words = name.split(" ")

        for word in words:
            # If adding this word would exceed the limit
            if current_line and len(current_line) + 1 + len(word) > 50:
                lines.append(current_line)
                current_line = word
            else:
                if current_line:
                    current_line += " " + word
                else:
                    current_line = word

        if current_line:
            lines.append(current_line)

        return "\n".join(lines)

    counted_micas_df['mica_name'] = counted_micas_df['mica_name'].map(mica_name_split)

    # convert counted_micas_df to heatmap format
    counted_micas_heatmap = pd.DataFrame(index=counted_micas_df['mica_name'], columns=disease_info['code_chapter'],
                                         dtype=int)

    for _, row in counted_micas_df.iterrows():
        counted_micas_heatmap.loc[row['mica_name'], row['diseases'].split('|')] = 1
    counted_micas_heatmap.fillna(0, inplace=True)

    top_n_pathways_counts = counted_micas_df.iloc[:n_top_pathways, :].copy()

    # Figure layout
    fig = plt.figure(figsize=(16, 12))
    fig.set_constrained_layout(True)
    gs = fig.add_gridspec(1, 3, width_ratios=[8, 1, 1], wspace=0.05)

    # --- Scatter plot ---
    ax = fig.add_subplot(gs[0])

    sns.heatmap(
        data=counted_micas_heatmap.iloc[:n_top_pathways, :], linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
        cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
        linecolor='black'
    )

    ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
    ax.set_ylabel('Reactome term')
    ax.set_xticklabels(ax.get_xticklabels(), size=10)
    ax.set_yticklabels(ax.get_yticklabels(), size=12)

    # --- Bar chart (counts per disease) ---
    # log_vars = np.log1p(n_vars_per_disease_full.values)
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    ax2.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['n_diseases'], color='grey',
             alpha=0.7)
    ax2.set_xlabel('n(diseases)')
    ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

    # Add vertical grid lines at tick marks
    xticks = ax2.get_xticks()
    for xt in xticks:
        ax2.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['n_diseases'].to_list()):
        ax2.text(count, y, str(count), va='center', fontsize=8)

    # --- Bar chart (gene counts per disease) ---
    # log_genes = np.log1p(n_genes_per_disease_full.values)
    ax3 = fig.add_subplot(gs[2], sharey=ax)
    ax3.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['pair_count'], color='grey',
             alpha=0.7)
    ax3.set_xlabel('n(disease pairs)')
    ax3.tick_params(axis='y', left=False, labelleft=False)

    # Vertical grid lines
    for xt in ax3.get_xticks():
        ax3.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['pair_count'].to_list()):
        ax3.text(count, y, str(count), va='center', fontsize=8)

    plt.show()

    plt.savefig(f"{plot_dir}/{plot_name}.png")
    plt.savefig(f"{plot_dir}/{plot_name}.svg")

    plt.close()

    return


""" Explore different filtering regimes """
# top 10% Resnik and contribution - biases towards low-level (high IC) terms
top_quantile_micas = bma_mica_contributions[
    (bma_mica_contributions['resnik'] >= bma_mica_contributions['resnik'].quantile(0.9)) &  # 90% = 8.202665
    (bma_mica_contributions['contribution'] >= bma_mica_contributions['contribution'].quantile(0.9))]  # 90% = 0.042463
top_quantile_micas_counted = count_high_flying_micas(top_quantile_micas)

# top 5 per disease pair by contribution (Resnik correlates with contribution quite well - sns.scatterplot(bma_mica_contributions, x='resnik', y='contribution'))
top5_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:5, :] for _, disease_pair_df in bma_mica_contributions.groupby(['disease1', 'disease2'])])
top5_micas_counted = count_high_flying_micas(top5_micas_per_disease_pair)

# top 5 per disease pair by contribution + MICA IC >= 25th percentile of all Reactome ICs
top5_thresholded_micas_per_disease_pair = top5_micas_per_disease_pair[
    top5_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
top5_thresholded_micas_counted = count_high_flying_micas(top5_thresholded_micas_per_disease_pair)

# identify complete connectors from the above filtering regime
top5_thresholded_micas_counted['pairwise_completion'] = [
    top5_thresholded_micas_counted.at[i, 'pair_count'] / math.comb(top5_thresholded_micas_counted.at[i, 'n_diseases'],
                                                                   2) for i in
    range(len(top5_thresholded_micas_counted))]
top5_thresholded_micas_complete_connectors = top5_thresholded_micas_counted[
    top5_thresholded_micas_counted['pairwise_completion'] == 1].drop(columns='pairwise_completion').copy()
top5_thresholded_micas_complete_connectors = top5_thresholded_micas_complete_connectors[
    top5_thresholded_micas_complete_connectors['n_diseases'] > 2].copy()

# plot
plot_dir = "~/ch2_genomics/2.12 pathway_based_similarities/plots"
plot_name = f"semsim_common_driver_pathways_reactome_absor_90perc_lifestyles"
plot_high_flying_micas(counted_micas_df=top_quantile_micas_counted, n_top_pathways=25, plot_name=plot_name, plot_dir=plot_dir)
plot_name = f"semsim_common_driver_pathways_reactome_absor_top5micas_lifestyles"
plot_high_flying_micas(counted_micas_df=top5_micas_counted, n_top_pathways=25, plot_name=plot_name, plot_dir=plot_dir)
plot_name = f"semsim_common_driver_pathways_reactome_absor_top5micas_thresholded_lifestyles"
plot_high_flying_micas(counted_micas_df=top5_thresholded_micas_counted, n_top_pathways=25, plot_name=plot_name, plot_dir=plot_dir)
plot_name = f"semsim_common_driver_pathways_reactome_absor_top5micas_thresholded_complete_connectors_lifestyles"
plot_high_flying_micas(counted_micas_df=top5_thresholded_micas_complete_connectors, n_top_pathways=25,
                       plot_name=plot_name, plot_dir=plot_dir)

""" Plot different IC/contribution dists for each filtering regime """
filtering = [
    'Unfiltered',
    'IC > 90th;\nContrib. > 90th',
    'Top-5',
    'Top-5; \nIC > global 25th'
]
dists = [bma_mica_contributions, top_quantile_micas, top5_micas_per_disease_pair, top5_thresholded_micas_per_disease_pair]
dists_to_plot = pd.concat([x[['resnik', 'contribution']].assign(filtering=filtering[i]) for i, x in enumerate(dists)])
dist_pairs = list(itertools.combinations(filtering, 2))
resnik_pvalues = [scipy.stats.mannwhitneyu(x[0].resnik, x[1].resnik).pvalue for x in list(itertools.combinations(dists, 2))]
contribution_pvalues = [scipy.stats.mannwhitneyu(x[0].contribution, x[1].contribution).pvalue for x in list(itertools.combinations(dists, 2))]

fig = plt.figure(figsize=(12, 8))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)

# --- Resnik (IC) ---
ax = fig.add_subplot(gs[0])
sns.violinplot(data=dists_to_plot, x='filtering', y='resnik', hue='filtering', hue_order=filtering, ax=ax, cut=0)
ax.axhline(0, color='black', linestyle='--', linewidth=0.5)
ax.set_xlabel('Filtering regime')
ax.set_ylabel('Information content')
# statistical annotations
annotator = Annotator(ax, dist_pairs, data=dists_to_plot, x='filtering', y='resnik')
annotator.configure(test_short_name='MW', text_format='full', loc='outside')
annotator.set_pvalues(resnik_pvalues)
annotator.annotate()
# ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.25)

# --- Contribution ---
ax1 = fig.add_subplot(gs[1])
sns.violinplot(data=dists_to_plot, x='filtering', y='contribution', hue='filtering', hue_order=filtering, ax=ax1, cut=0)
ax1.axhline(0, color='black', linestyle='--', linewidth=0.5)
ax1.set_xlabel('Filtering regime')
ax1.set_ylabel('Contribution')
# statistical annotations
annotator = Annotator(ax1, dist_pairs, data=dists_to_plot, x='filtering', y='contribution')
annotator.configure(test_short_name='MW', text_format='full', loc='outside')
annotator.set_pvalues(contribution_pvalues)
annotator.annotate()

# export
ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.6)
ax1.set_ylim(ax1.get_ylim()[0], ax1.get_ylim()[1] * 1.6)
plt.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/semsim_driver_pathways_reactome_dists.svg")
plt.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/semsim_driver_pathways_reactome_dists.png")
plt.close()

""" Identifying driver pathways with the final filtering regime """
# top 5 per disease pair by contribution (Resnik correlates with contribution quite well - sns.scatterplot(bma_mica_contributions, x='resnik', y='contribution'))
top5_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:5, :] for _, disease_pair_df in bma_mica_contributions.groupby(['disease1', 'disease2'])])
top5_micas_counted = count_high_flying_micas(top5_micas_per_disease_pair)

# top 5 per disease pair by contribution + MICA IC >= 25th percentile of all Reactome ICs
top5_thresholded_micas_per_disease_pair = top5_micas_per_disease_pair[
    top5_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
top5_thresholded_micas_counted = count_high_flying_micas(top5_thresholded_micas_per_disease_pair)

# identify complete connectors
top5_thresholded_micas_counted['pairwise_completion'] = [
    top5_thresholded_micas_counted.at[i, 'pair_count'] / math.comb(top5_thresholded_micas_counted.at[i, 'n_diseases'],
                                                                   2) for i in
    range(len(top5_thresholded_micas_counted))]
top5_thresholded_micas_complete_connectors = top5_thresholded_micas_counted[
    top5_thresholded_micas_counted['pairwise_completion'] == 1].drop(columns='pairwise_completion').copy()
top5_thresholded_micas_complete_connectors = top5_thresholded_micas_complete_connectors[
    top5_thresholded_micas_complete_connectors['n_diseases'] > 2].copy()
top5_thresholded_micas_counted.to_csv(
    f"~/data/internal/genomics/pathway_based_similarities/semsim_common_driver_pathways_reactome_absor_top5micas_thresholded.csv",
    index=False)
top5_thresholded_micas_complete_connectors.to_csv(
    f"~/data/internal/genomics/pathway_based_similarities/semsim_common_driver_pathways_reactome_absor_top5micas_thresholded_complete_connectors.csv",
    index=False)

# parentage analysis
from utils.enrichment_analyses import reactome_parentage_analysis
_ = reactome_parentage_analysis(annotation_set=top5_thresholded_micas_complete_connectors.mica.to_list(),
                                plot_path=f"~/ch2_genomics/2.12 pathway_based_similarities/plots/semsim_common_driver_pathways_reactome_absor_top5micas_thresholded_complete_connectors_lifestyles_parentage",
                                export_plot=True)

""" FCA to get disease-driver MICA groups from driver MICA-based groups of diseases """
# Build incidence matrix: rows = groups, columns = diseases
# diseases = sorted({d for g in groups for d in g})
# df = pd.DataFrame(False, index=range(len(groups)), columns=diseases)
# for i, g in enumerate(groups):
#     df.loc[i, list(g)] = True
fca_prep_df = top5_thresholded_micas_complete_connectors[['mica', 'diseases']].copy()
fca_prep_df['diseases'] = fca_prep_df['diseases'].str.split('|')
fca_prep_df = fca_prep_df.explode('diseases').assign(values=True)
fca_mat_df = fca_prep_df.pivot(index='mica', columns='diseases', values='values').fillna(False)

# Run FCA
ctx = Context(objects=fca_mat_df.index.astype(str).tolist(),
              properties=fca_mat_df.columns.astype(str).tolist(),
              bools=fca_mat_df.astype(bool).values.tolist())
lattice = ctx.lattice

# Extract interpretable overlap patterns
overlap_patterns_list = []
bar = ProgBar(len(lattice), stream=sys.stdout, title='FCA for disease groups by pathways')
for concept in lattice:
    diseases_in_common = concept.intent  # exact disease set
    groups_sharing = concept.extent  # groups that share them

    # Filter trivial concepts
    if len(diseases_in_common) >= 3 and len(groups_sharing) >= 2:
        overlap_patterns_list.append({
            'n_diseases': len(diseases_in_common),
            'n_pathways': len(groups_sharing),
            'diseases': '|'.join(sorted(list(diseases_in_common))),
            'pathways': '|'.join(sorted(list(groups_sharing)))
        })
        # print(f"Diseases {set(diseases_in_common)}"
        #       f"appear in groups {set(groups_sharing)}")

    bar.update()
overlap_patterns_df = pd.DataFrame(overlap_patterns_list)
overlap_patterns_df['rank_measure'] = overlap_patterns_df['n_diseases'] * overlap_patterns_df['n_pathways']
overlap_patterns_df = overlap_patterns_df.sort_values(by='rank_measure', ascending=False)
overlap_patterns_df.to_csv(f"~/data/internal/genomics/pathway_based_similarities/fca_results.csv", index=False)

# add Reactome names to FCA
# get reactome pathway names
reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                header=None, names=['id', 'name', 'species'])
reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
id_to_name_map = dict(zip(reactome_pathways.id, reactome_pathways.name))
overlap_patterns_for_export = overlap_patterns_df.drop(columns='rank_measure')
overlap_patterns_for_export['pathways'] = overlap_patterns_for_export['pathways'].apply(lambda x: '; '.join([id_to_name_map[pathway] for pathway in x.split('|')]))
overlap_patterns_for_export.to_csv(f"~/data/internal/genomics/pathway_based_similarities/fca_results_named.csv", index=False)

pathway_lst = []
pathway_lst_trios = []
for overlap_pattern in overlap_patterns_list:
    if overlap_pattern['n_diseases'] == 3:
        pathway_lst_trios += overlap_pattern['pathways'].split('|')
    pathway_lst += overlap_pattern['pathways'].split('|')
overlap_pathway_freq = pd.Series(pathway_lst).value_counts()
overlap_pathway_trios_freq = pd.Series(pathway_lst_trios).value_counts()
len(set(pathway_lst))

# plot pathway frequency in all sets and only in trio sets
freqs_to_plot = pd.concat([pd.DataFrame(overlap_pathway_freq).reset_index().set_axis(['pathway', 'n(sets)'], axis=1).assign(set_type='All sets'),
                           pd.DataFrame(overlap_pathway_trios_freq).reset_index().set_axis(['pathway', 'n(sets)'], axis=1).assign(set_type='Trios only')])
freqs_to_plot = freqs_to_plot.merge(reactome_pathways.set_axis(['pathway', 'Reactome term'], axis=1))

g = sns.catplot(data=freqs_to_plot, x='n(sets)', y='Reactome term', col='set_type', kind='bar', sharey=True,
                sharex=True, height=6, aspect=0.67)
g.set_axis_labels('n(sets)', 'Reactome term')
g.set_titles("{col_name}", size=10)
g.axes.flat[0].tick_params(axis='y', labelsize=5)
for ax in g.axes.flat:
    ax.tick_params(axis='x', labelsize=8)
    ax.set_axisbelow(True)   # bands/lines go behind bars
    yticks = ax.get_yticks()
    # alternate shading
    for i in range(len(yticks)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, alpha=0.1, color='lightgrey')
    # add horizontal lines between categories
    for i in range(len(yticks) - 1):
        ax.axhline(i + 0.5, color='grey', linewidth=0.5, alpha=0.5)
    ax.set_ylim(-0.6, len(yticks) - 0.4)  # remove additional padding, while retaining a thin one
    ax.invert_yaxis()  # above line inverts y-axis, so un-invert
    ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
g.tight_layout()
g.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/fca_reactome_disease_set_trio_freqs.svg")
g.savefig(f"~/ch2_genomics/2.12 pathway_based_similarities/plots/fca_reactome_disease_set_trio_freqs.png")
plt.close()
