"""
Analyse results from embedding similarity tests to explore dislocalisation between GWAS and proteomics genes.
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
import numpy as np

from utils.significance_labelling import significance_labelling

mpl.use('TkAgg')

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Read in all results data """
results_list = []
for dirr in [f"~/data/internal/knowledge_graph/gene_protein_distances/embeddings/{x.name}" for x in os.scandir('~/data/internal/knowledge_graph/gene_protein_distances/embeddings/') if x.is_dir()]:
    results_list = results_list + [pd.read_csv(f"{dirr}/{x.name}").assign(proteomics=x.name.split('.')[0].split('_')[-1],
                                                                          gene_gene_links='_'.join(x.name.split('_')[2:-1])) for x in os.scandir(dirr) if x.is_file()]
results = pd.concat(results_list)

""" Summary plot """
results_long = results[[x for x in results.columns if 'biserial' not in x]].melt(id_vars=['disease', 'tissue', 'gwas_gene_disease_relationship', 'proteomics', 'gene_gene_links'],
                                                                                 value_vars=[x for x in results.columns if 'inter_pval' in x], value_name='pval', var_name='comparison')
results_long['comparison'] = results_long['comparison'].replace({'intra_gwas_vs_inter_pval': 'Intra-GWAS vs Inter',
                                                             'intra_proteomics_vs_inter_pval': 'Intra-proteomics vs Inter'})
results_long_select = results_long[results_long['pval'].notna()]
results_long_select = significance_labelling(results_long_select, pvalue_col='pval')
results_long_select['x_axis_category'] = results_long_select.proteomics + '|' + results_long_select.gene_gene_links + '|' + results_long_select.gwas_gene_disease_relationship

len_each_x_axis_cat = pd.concat([pd.DataFrame(sub_df.comparison.value_counts()).reset_index().assign(x_axis_category=cat)
                             for cat, sub_df in results_long_select.groupby('x_axis_category')])

results_long_select_sig = results_long_select[results_long_select.fdr_sig == 1].copy()

summary_plot_df = pd.concat([pd.DataFrame(sub_df.comparison.value_counts()).reset_index().assign(x_axis_category=cat)
                             for cat, sub_df in results_long_select_sig.groupby('x_axis_category')])
summary_plot_df = summary_plot_df.merge(len_each_x_axis_cat.rename(columns={'count': 'total_count'}))
summary_plot_df['proportion'] = summary_plot_df['count'] * 100 / summary_plot_df['total_count']

# stats for text
summary_plot_df[summary_plot_df['comparison'].str.contains('proteomics')].proportion.describe()
summary_plot_df[summary_plot_df['comparison'].str.contains('GWAS')].proportion.describe()

plt.figure(figsize=(12, 8))
sns.barplot(data=summary_plot_df.rename(columns={'comparison': 'Comparison'}), x='x_axis_category', y='proportion',
            hue='Comparison',
            hue_order=['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter'])
plt.xlabel('Options used in graph building (GWAS-gene links | proteomics algorithm | gene-gene links)')
plt.ylabel('Proportion of tissue-ARD pairs')
plt.gca().tick_params(axis='x', rotation=90)
plt.grid(True, which='major', axis='y', color='lightgrey', lw=0.4, ls='--')
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig('~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/overall_embeddings_summary.png')
plt.savefig('~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/overall_embeddings_summary.svg')
plt.close()

""" Scattermaps with x/y-axis barplots for each graph-build option """
# map disease data to categorical index space
all_disease_labels = sorted(set(results.merge(disease_info[['icd10_three_letter', 'code_chapter']].rename(columns={'icd10_three_letter': 'disease'})).code_chapter))
# make mapping have 1-tick gaps between chapters, so that there are whitespace gaps/boundaries between chapters
disease_label_to_num = {}
current_pos = 0
for i in range(len(all_disease_labels)):
    if i == 0:
        disease_label_to_num[all_disease_labels[i]] = current_pos
        current_pos += 1
        continue
    current_chapter = int(re.search(r"\((\d+)\)", all_disease_labels[i]).group(1))
    past_chapter = int(re.search(r"\((\d+)\)", all_disease_labels[i - 1]).group(1))
    if past_chapter < current_chapter:
        current_pos += 1  # insert a 1-tick wide gap into the axis at the boundary between chapters
    disease_label_to_num[all_disease_labels[i]] = current_pos
    current_pos += 1
# label_to_num = {lab: i for i, lab in enumerate(all_labels)}
code_chapter_to_code = dict(zip(disease_info.code_chapter, disease_info.icd10_three_letter))
disease_code_to_num = {code_chapter_to_code[code_chapter]: num for code_chapter, num in disease_label_to_num.items()}

# map tissue data to categorical index space
all_tissue_labels = sorted(set(results.tissue))
tissue_label_to_num = {e: i for i, e in enumerate(all_tissue_labels)}

results_long_select_sig_slim = results_long_select_sig[['disease', 'tissue', 'comparison', 'x_axis_category']].copy()
results_long_select_sig_slim['disease_as_num'] = results_long_select_sig_slim['disease'].map(disease_code_to_num)
results_long_select_sig_slim['tissue_as_num'] = results_long_select_sig_slim['tissue'].map(tissue_label_to_num)

mpl.use('Agg')
for graph_build_options in results_long_select_sig_slim.x_axis_category.unique():
    # for graph_build_options in ...:
    # single_result = results_long_select_sig_slim[results_long_select_sig_slim.x_axis_category == 'magma-linked|firth|coexpresses-with']
    single_result = results_long_select_sig_slim[results_long_select_sig_slim.x_axis_category == graph_build_options]
    results_pie_points = [{'tissue_as_num': tissue_disease_tup[0], 'disease_as_num': tissue_disease_tup[1],
                           'fracs': list(tissue_disease_df.comparison.value_counts().reindex(
                               index=['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']).fillna(0) / len(tissue_disease_df))}
                          for tissue_disease_tup, tissue_disease_df in single_result.groupby(['tissue_as_num', 'disease_as_num'])]

    # --- Sample data ---
    # Each point has an (x, y) position and a composition (list of fractions)
    # points = [
    #     {"x": 1, "y": 2, "fracs": [0.5, 0.3, 0.2], "label": "A"},
    #     {"x": 3, "y": 5, "fracs": [0.2, 0.5, 0.3], "label": "B"},
    #     {"x": 5, "y": 1, "fracs": [0.1, 0.1, 0.8], "label": "C"},
    #     {"x": 7, "y": 4, "fracs": [0.4, 0.4, 0.2], "label": "D"},
    #     {"x": 2, "y": 7, "fracs": [0.6, 0.2, 0.2], "label": "E"},
    # ]

    # 1st 2 tab10 colours
    colors = ["#1f77b4", "#ff7f0e"]
    colors_dict = {'Intra-GWAS vs Inter': '#1f77b4', 'Intra-proteomics vs Inter': '#ff7f0e'}
    pie_size = 0.01  # size of each pie in figure-fraction units

    fig = plt.figure(figsize=(12, 16), dpi=100)
    fig.set_constrained_layout(True)
    fig.set_constrained_layout_pads(w_pad=0.1)
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 6], wspace=0.05)

    # --- pie chart scatter ---
    ax2 = fig.add_subplot(gs[2])
    # Set axis limits with padding so pies aren't clipped at edges
    xs = [p['tissue_as_num'] for p in results_pie_points]
    ys = [p['disease_as_num'] for p in results_pie_points]
    ax2.set_xlim(min(xs) - 1.5, max(xs) + 1.5)
    ax2.set_ylim(min(list(disease_code_to_num.values())) - 1.5, max(list(disease_code_to_num.values())) + 1.5)
    ax2.set_xlabel('Tissue')
    ax2.set_ylabel('Disease ICD-10 code (ICD-10 chapter)')

    # if ax2.get_yaxis().get_inverted():
    ax2.invert_yaxis()

    # set ticks labels to the string versions
    x_tick_positions = [tissue_label_to_num[label] for label in all_tissue_labels]
    ax2.set_xticks(x_tick_positions)
    ax2.set_xticklabels(all_tissue_labels, rotation=90, size=7)
    y_tick_positions = [disease_label_to_num[label] for label in all_disease_labels]
    ax2.set_yticks(y_tick_positions)
    ax2.set_yticklabels(all_disease_labels, size=7)
    plt.grid(True, zorder=0, color='lightgrey', lw=0.4, ls='--')

    inset_axes_list = []  # keep references so we can reposition them


    def place_pies(event=None):
        # Remove previously placed inset axes (on redraw)
        for ia in inset_axes_list:
            ia.remove()
        inset_axes_list.clear()

        for point in results_pie_points:
            # Convert data coords → figure-fraction coords
            disp = ax2.transData.transform((point['tissue_as_num'], point['disease_as_num']))
            fig_coord = fig.transFigure.inverted().transform(disp)

            # Place an inset axis centred on the point
            inset_ax = fig.add_axes(
                [
                    fig_coord[0] - pie_size / 2,
                    fig_coord[1] - pie_size / 2,
                    pie_size,
                    pie_size,
                ]
            )
            inset_ax.pie(point["fracs"], colors=colors)
            inset_ax.set_aspect("equal")
            inset_axes_list.append(inset_ax)


    # --- number of non-neither categories per tissue, over all diseases ---
    n_cats_per_tissue = pd.DataFrame([{'tissue': tissue} | tissue_df.comparison.value_counts().to_dict()
                                      for tissue, tissue_df in single_result.groupby('tissue')])
    for cat in ['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']:
        if cat not in n_cats_per_tissue.columns:
            n_cats_per_tissue[cat] = 0
    n_cats_per_tissue_long = n_cats_per_tissue.melt(id_vars='tissue', var_name='category', value_name='counts')
    n_cats_per_tissue_long['log_count'] = np.log1p(n_cats_per_tissue_long['counts'])
    n_cats_per_tissue_long['tissue_as_num'] = n_cats_per_tissue_long['tissue'].map(tissue_label_to_num)

    ax0 = fig.add_subplot(gs[0], sharex=ax2)
    bar_width = 0.25
    offsets = {'Intra-GWAS vs Inter': -bar_width/2, 'Intra-proteomics vs Inter': bar_width/2}
    for cat in ['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']:
        cat_data = n_cats_per_tissue_long[n_cats_per_tissue_long['category'] == cat]
        ax0.bar(
            cat_data['tissue_as_num'] + offsets[cat],
            cat_data['log_count'],
            width=bar_width,
            color=colors_dict[cat],
            label=cat
        )
    ax0.set_xlim(ax2.get_xlim())
    ax0.set_ylabel('ln(count + 1)')
    ax0.set_xlabel('')
    ax0.tick_params(axis='x', which='both', bottom=False, labelbottom=False)  # hide y-axis labels on RHS plot

    # Add horizontal grid lines at tick marks
    ax0.grid(True, which='major', axis='y', color='lightgrey', lw=0.4, ls='--')
    ax0.set_axisbelow(True)

    # --- number of non-neither categories per disease, over all tissues ---
    n_cats_per_disease = pd.DataFrame([{'disease': disease} | disease_df.comparison.value_counts().to_dict()
                                       for disease, disease_df in
                                       single_result.groupby('disease')]).fillna(0)
    for cat in ['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']:
        if cat not in n_cats_per_disease.columns:
            n_cats_per_disease[cat] = 0
    n_cats_per_disease_long = n_cats_per_disease.melt(id_vars='disease', var_name='category', value_name='counts')
    n_cats_per_disease_long['log_count'] = np.log1p(n_cats_per_disease_long['counts'])
    n_cats_per_disease_long['disease_as_num'] = n_cats_per_disease_long['disease'].map(disease_code_to_num)

    ax3 = fig.add_subplot(gs[3])
    bar_height = 0.25
    for cat in ['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']:
        cat_data = n_cats_per_disease_long[n_cats_per_disease_long['category'] == cat]
        ax3.barh(
            cat_data['disease_as_num'] + offsets[cat],
            cat_data['log_count'],
            height=bar_height,
            color=colors_dict[cat],
            label=cat
        )

    # sns.barplot(data=n_cats_per_disease_long, y='disease_as_num', x='log_count', hue='category', palette='tab10',
    #             hue_order=['both', 'gwas', 'proteomics'], ax=ax3)
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_xlabel('ln(count + 1)')
    ax3.set_ylabel('')
    ax3.tick_params(axis='y', which='both', left=False, labelleft=False)  # hide y-axis labels on RHS plot

    # Add horizontal grid lines at tick marks
    ax3.grid(True, which='major', axis='x', color='lightgrey', lw=0.4, ls='--')
    ax3.set_axisbelow(True)

    # --- draw legend in top-RH corner (gs[1])
    ax1 = fig.add_subplot(gs[1])
    legend_labels = ['Intra-GWAS vs Inter', 'Intra-proteomics vs Inter']
    ax1.legend(handles=[mpl.patches.Patch(color=color, label=label) for label, color in
                        dict(zip(legend_labels, colors)).items()],
               loc='center')
    ax1.axis('off')

    # --- draw pies ---
    # Fire once after the first layout pass, then on every subsequent redraw
    fig.canvas.mpl_connect("draw_event", place_pies)
    # place_pies()
    # plt.pause(0.1)
    fig.canvas.draw()
    # plt.show()
    # fig.canvas.draw()

    fig.savefig(
        f"~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/n_sig_comparisons_per_combo_{graph_build_options.replace('|', '_')}.png",
        dpi=100)
    fig.savefig(
        f"~/ch5_integration/5.4 graph_localisation/5.4.2 embeddings/plots/n_sig_comparisons_per_combo_{graph_build_options.replace('|', '_')}.svg",
        dpi=100)

    plt.close()
