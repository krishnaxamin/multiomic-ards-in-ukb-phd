"""
Plot cluster enrichment results to show
(1) the number of clusters across each tissue-disease combination enriched in either GWAS, proteomics or both
(2) enrichment for co-occurrence of GWAS&proteomics enrichment ('both' from (1)) in each tissue-disease combination
"""
import pandas as pd
import os
import re
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')

disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Read in results data """
# result files contain per-cluster GWAS-genes and proteomics-genes enrichment for each tissue-disease combination
# enrichment stats and significance labels for each of GWAS-genes and proteomics-genes are in 'gwas_*' and
#  'proteomics_*' respectively
# specificity labels at FDR-sig or nominal = whether cluster is enriched for one of, both, or neither of GWAS-genes and
#  proteomics-genes
results_files = [x.name for x in os.scandir(f"~/data/internal/knowledge_graph/gene_protein_distances/clustering")
                 if
                 x.is_file() and 'results' in x.name]
# summary files contain per-tissue-disease-combination enrichment in clusters that host both GWAS-genes and
#  proteomics-genes enrichment
# the per-tissue-disease-combination enrichment is performed for cluster-wise nominal and FDR-corrected significance
#  (col='sig_level')
# enrichment pvalues are in 'fisher_less' and 'fisher_greater', with significance labels in '*nom_sig' and '*fdr_sig'
summaries_files = [x.name for x in
                   os.scandir(f"~/data/internal/knowledge_graph/gene_protein_distances/clustering") if
                   x.is_file() and 'summaries' in x.name]

results_dict = {re.search(r'(?<=results_)[^.]+', result_file_name).group():
    pd.read_csv(
        f"~/data/internal/knowledge_graph/gene_protein_distances/clustering/{result_file_name}")
    for result_file_name in results_files}
summaries_dict = {re.search(r'(?<=summaries_)[^.]+', summary_file_name).group():
    pd.read_csv(
        f"~/data/internal/knowledge_graph/gene_protein_distances/clustering/{summary_file_name}")
    for summary_file_name in summaries_files}

assert list(results_dict.keys()) == list(
    summaries_dict.keys()), 'Difference sets of results in "results" and "summaries" file sets.'
analysis_categories = list(results_dict.keys())

""" Plot """

for analysis_category in analysis_categories:

    for gene_gwas_relationship in ['variant-linked', 'magma-linked', 'variant-linked_magma-linked']:
        analysis_category_with_genegwas = f"{analysis_category}_{gene_gwas_relationship}"

        print(analysis_category_with_genegwas)
        algo = analysis_category.split('_')[-1]

        """ Plot proportion of clusters FDR-enriched in GWAS-genes, proteomics-genes, or both, per tissue-disease """
        print('- graph1')
        # analysis_category = 'coexpresses-with_cox'
        cluster_prop_df = results_dict[analysis_category].copy()
        cluster_prop_df = cluster_prop_df[cluster_prop_df['gwas_gene_disease_relationship'] == gene_gwas_relationship]
        cluster_prop_df = (cluster_prop_df.merge(disease_info[['icd10_three_letter', 'code_chapter']]
                                                 .rename(columns={'icd10_three_letter': 'disease'}))
                           .drop(columns='disease').rename(columns={'code_chapter': 'disease'}))
        cluster_prop_df_no_neither = cluster_prop_df[cluster_prop_df.fdr_specificity_label != 'neither'].copy()

        # map disease data to categorical index space
        all_disease_labels = sorted(set(cluster_prop_df.disease))
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
        cluster_prop_df_no_neither['disease_as_num'] = cluster_prop_df_no_neither['disease'].map(disease_label_to_num)

        # map tissue data to categorical index space
        all_tissue_labels = sorted(set(cluster_prop_df.tissue))
        tissue_label_to_num = {e: i for i, e in enumerate(all_tissue_labels)}
        cluster_prop_df_no_neither['tissue_as_num'] = cluster_prop_df_no_neither['tissue'].map(tissue_label_to_num)

        # fracs: [both, gwas-only, proteomics-only, (neither)]
        # fracs: if 'neither' not included, then the pie chart covers only those tissue-disease combos with enrichment
        cluster_prop_df_points = [{'tissue_as_num': tissue_disease_tup[0], 'disease_as_num': tissue_disease_tup[1],
                                   'fracs': list(tissue_disease_df.fdr_specificity_label.value_counts().reindex(
                                       index=['both', 'gwas', 'proteomics']).fillna(0) / len(tissue_disease_df))}
                                  for tissue_disease_tup, tissue_disease_df in
                                  cluster_prop_df_no_neither.groupby(['tissue_as_num', 'disease_as_num'])]

        # --- Sample data ---
        # Each point has an (x, y) position and a composition (list of fractions)
        # points = [
        #     {"x": 1, "y": 2, "fracs": [0.5, 0.3, 0.2], "label": "A"},
        #     {"x": 3, "y": 5, "fracs": [0.2, 0.5, 0.3], "label": "B"},
        #     {"x": 5, "y": 1, "fracs": [0.1, 0.1, 0.8], "label": "C"},
        #     {"x": 7, "y": 4, "fracs": [0.4, 0.4, 0.2], "label": "D"},
        #     {"x": 2, "y": 7, "fracs": [0.6, 0.2, 0.2], "label": "E"},
        # ]

        # 1st three tab10 colours
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
        colors_dict = {'both': '#1f77b4', 'gwas': '#ff7f0e', 'proteomics': '#2ca02c'}
        pie_size = 0.01  # size of each pie in figure-fraction units

        fig = plt.figure(figsize=(12, 16), dpi=100)
        fig.set_constrained_layout(True)
        fig.set_constrained_layout_pads(w_pad=0.1)
        gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 6], wspace=0.05)

        # --- pie chart scatter ---
        ax2 = fig.add_subplot(gs[2])
        # Set axis limits with padding so pies aren't clipped at edges
        xs = [p['tissue_as_num'] for p in cluster_prop_df_points]
        ys = [p['disease_as_num'] for p in cluster_prop_df_points]
        ax2.set_xlim(min(xs) - 1.5, max(xs) + 1.5)
        ax2.set_ylim(min(ys) - 1.5, max(ys) + 1.5)
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

            for point in cluster_prop_df_points:
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
        n_cats_per_tissue = pd.DataFrame([{'tissue': tissue} | tissue_df.fdr_specificity_label.value_counts().to_dict()
                                          for tissue, tissue_df in cluster_prop_df_no_neither.groupby('tissue')])
        for cat in ['proteomics', 'gwas', 'both']:
            if cat not in n_cats_per_tissue.columns:
                n_cats_per_tissue[cat] = 0
        n_cats_per_tissue_long = n_cats_per_tissue.melt(id_vars='tissue', var_name='category', value_name='counts')
        n_cats_per_tissue_long['log_count'] = np.log1p(n_cats_per_tissue_long['counts'])
        n_cats_per_tissue_long['tissue_as_num'] = n_cats_per_tissue_long['tissue'].map(tissue_label_to_num)

        ax0 = fig.add_subplot(gs[0], sharex=ax2)
        bar_width = 0.25
        offsets = {'both': -bar_width, 'gwas': 0, 'proteomics': bar_width}
        for cat in ['both', 'gwas', 'proteomics']:
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
        n_cats_per_disease = pd.DataFrame(
            [{'disease': disease} | disease_df.fdr_specificity_label.value_counts().to_dict()
             for disease, disease_df in
             cluster_prop_df_no_neither.groupby('disease')]).fillna(0)
        for cat in ['proteomics', 'gwas', 'both']:
            if cat not in n_cats_per_disease.columns:
                n_cats_per_disease[cat] = 0
        n_cats_per_disease_long = n_cats_per_disease.melt(id_vars='disease', var_name='category', value_name='counts')
        n_cats_per_disease_long['log_count'] = np.log1p(n_cats_per_disease_long['counts'])
        n_cats_per_disease_long['disease_as_num'] = n_cats_per_disease_long['disease'].map(disease_label_to_num)

        ax3 = fig.add_subplot(gs[3])
        bar_height = 0.25
        for cat in ['both', 'gwas', 'proteomics']:
            cat_data = n_cats_per_disease_long[n_cats_per_disease_long['category'] == cat]
            ax3.barh(
                cat_data['disease_as_num'] + offsets[cat],
                cat_data['log_count'],
                height=bar_height,
                color=colors_dict[cat],
                label=cat
            )

        ax3.set_ylim(ax2.get_ylim())
        ax3.set_xlabel('ln(count + 1)')
        ax3.set_ylabel('')
        ax3.tick_params(axis='y', which='both', left=False, labelleft=False)  # hide y-axis labels on RHS plot

        # Add horizontal grid lines at tick marks
        ax3.grid(True, which='major', axis='x', color='lightgrey', lw=0.4, ls='--')
        ax3.set_axisbelow(True)

        # --- draw legend in top-RH corner (gs[1])
        ax1 = fig.add_subplot(gs[1])
        legend_labels = ["Both", "GWAS only", f"Proteomics ({algo}) only"]
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
            f"~/ch5_integration/5.4 graph_localisation/plots/n_clusters_per_combo_{analysis_category_with_genegwas}.png",
            dpi=100)
        fig.savefig(
            f"~/ch5_integration/5.4 graph_localisation/plots/n_clusters_per_combo_{analysis_category_with_genegwas}.svg",
            dpi=100)

        plt.close()

        """ Plot enrichment of co-occurrence per tissue-disease combination """
        print('- graph2')
        # analysis_category = 'coexpresses-with_cox'
        cooccurence_df = summaries_dict[analysis_category].copy()
        cooccurence_df = cooccurence_df[cooccurence_df['gwas_gene_disease_relationship'] == gene_gwas_relationship]
        cooccurence_df = (cooccurence_df.merge(disease_info[['icd10_three_letter', 'code_chapter']]
                                               .rename(columns={'icd10_three_letter': 'disease'}))
                          .drop(columns='disease').rename(columns={'code_chapter': 'disease'}))
        cooccurence_df = cooccurence_df[cooccurence_df.sig_level == 'fdr'][
            ['tissue', 'disease', 'fisher_less_nom_sig', 'fisher_greater_nom_sig',
             'fisher_less_fdr_sig', 'fisher_greater_fdr_sig']].copy()

        cooccurence_df_long = cooccurence_df.melt(id_vars=['tissue', 'disease'], value_name='counts',
                                                  var_name='category')

        # map disease data to categorical index space
        all_disease_labels = sorted(set(cooccurence_df_long.disease))
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
        cooccurence_df_long['disease_as_num'] = cooccurence_df_long['disease'].map(disease_label_to_num)

        # map tissue data to categorical index space
        all_tissue_labels = sorted(set(cooccurence_df_long.tissue))
        tissue_label_to_num = {e: i for i, e in enumerate(all_tissue_labels)}
        cooccurence_df_long['tissue_as_num'] = cooccurence_df_long['tissue'].map(tissue_label_to_num)

        cooccurence_df_long = cooccurence_df_long[cooccurence_df_long.counts > 0].copy()
        if len(cooccurence_df_long) == 0:
            print(
                f"graph 2 for {analysis_category}_{gene_gwas_relationship} does not have significant enrichment for co-occurrence at any level.")
            continue
        # keep the most significant entry for each tissue-disease combination
        cooccurence_df_long_list = []
        for tissue_disease_tup, tissue_disease_df in cooccurence_df_long.groupby(['tissue_as_num', 'disease_as_num']):
            if len(tissue_disease_df) > 1:
                cooccurence_df_long_list.append(
                    tissue_disease_df[tissue_disease_df.category.str.contains('fdr_sig')].copy())
            else:
                cooccurence_df_long_list.append(tissue_disease_df)
        cooccurence_df_long_most_sig = pd.concat(cooccurence_df_long_list)

        # fdr-sig greater: plt.coolwarm(1.0). fdr-sig less: plt.coolwarm(0.0)
        # nom-sig greater: paler version of plt.coolwarm(1.0). nom-sig less: paler version of plt.coolwarm(0.0)
        colors_dict = {'fisher_greater_fdr_sig': plt.cm.coolwarm(1.0),
                       'fisher_greater_nom_sig': '#FEB9C7',
                       'fisher_less_fdr_sig': plt.cm.coolwarm(0.0),
                       'fisher_less_nom_sig': '#C8CDEE'}

        fig = plt.figure(figsize=(12, 16), dpi=100)
        fig.set_constrained_layout(True)
        fig.set_constrained_layout_pads(w_pad=0.1)
        gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 6], wspace=0.05)

        # --- scatter ---
        ax2 = fig.add_subplot(gs[2])
        sns.scatterplot(data=cooccurence_df_long_most_sig, x='tissue_as_num', y='disease_as_num',
                        hue='category', palette=colors_dict, legend=False)
        ax2.set_xlim(-0.55, max(list(tissue_label_to_num.values())) + 0.5)
        ax2.set_ylim(-0.5, max(list(disease_label_to_num.values())) + 0.5)
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
        ax2.grid(True, zorder=0, color='lightgrey', lw=0.4, ls='--')
        ax2.set_axisbelow(True)

        # --- number of non-neither categories per tissue, over all diseases ---
        n_cats_per_tissue = pd.DataFrame([{'tissue': tissue} | tissue_df.category.value_counts().to_dict()
                                          for tissue, tissue_df in cooccurence_df_long_most_sig.groupby('tissue')])
        for cat in colors_dict.keys():
            if cat not in n_cats_per_tissue.columns:
                n_cats_per_tissue[cat] = 0
        n_cats_per_tissue_long = n_cats_per_tissue.melt(id_vars='tissue', var_name='category', value_name='counts')
        n_cats_per_tissue_long['log_count'] = np.log1p(n_cats_per_tissue_long['counts'])
        n_cats_per_tissue_long['tissue_as_num'] = n_cats_per_tissue_long['tissue'].map(tissue_label_to_num)

        ax0 = fig.add_subplot(gs[0], sharex=ax2)
        bar_height_width = 0.15
        offsets = {'fisher_greater_fdr_sig': -bar_height_width * 1.5,
                   'fisher_greater_nom_sig': -bar_height_width * 0.5,
                   'fisher_less_fdr_sig': bar_height_width * 0.5,
                   'fisher_less_nom_sig': bar_height_width * 1.5}
        for cat in offsets.keys():
            cat_data = n_cats_per_tissue_long[n_cats_per_tissue_long['category'] == cat]
            ax0.bar(
                cat_data['tissue_as_num'] + offsets[cat],
                cat_data['log_count'],
                width=bar_height_width,
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
        n_cats_per_disease = pd.DataFrame([{'disease': disease} | disease_df.category.value_counts().to_dict()
                                           for disease, disease_df in
                                           cooccurence_df_long_most_sig.groupby('disease')]).fillna(0)
        for cat in colors_dict.keys():
            if cat not in n_cats_per_disease.columns:
                n_cats_per_disease[cat] = 0
        n_cats_per_disease_long = n_cats_per_disease.melt(id_vars='disease', var_name='category', value_name='counts')
        n_cats_per_disease_long['log_count'] = np.log1p(n_cats_per_disease_long['counts'])
        n_cats_per_disease_long['disease_as_num'] = n_cats_per_disease_long['disease'].map(disease_label_to_num)

        ax3 = fig.add_subplot(gs[3])
        bar_height = 0.25
        for cat in offsets.keys():
            cat_data = n_cats_per_disease_long[n_cats_per_disease_long['category'] == cat]
            ax3.barh(
                cat_data['disease_as_num'] + offsets[cat],
                cat_data['log_count'],
                height=bar_height_width,
                color=colors_dict[cat],
                label=cat
            )

        ax3.set_ylim(ax2.get_ylim())
        ax3.set_xlabel('ln(count + 1)')
        ax3.set_ylabel('')
        ax3.tick_params(axis='y', which='both', left=False, labelleft=False)  # hide y-axis labels on RHS plot

        # Add horizontal grid lines at tick marks
        ax3.grid(True, which='major', axis='x', color='lightgrey', lw=0.4, ls='--')
        ax3.set_axisbelow(True)

        # --- draw legend in top-RH corner (gs[1])
        ax1 = fig.add_subplot(gs[1])
        ax1.legend(handles=[
            mpl.patches.Patch(color=plt.cm.coolwarm(1.0), label='Greater (FDR < 0.05)'),
            mpl.patches.Patch(color='#FEB9C7', label='Greater (p < 0.05)'),
            mpl.patches.Patch(color=plt.cm.coolwarm(0.0), label='Less (FDR < 0.05)'),
            mpl.patches.Patch(color='#C8CDEE', label='Less (p < 0.05)'),
        ], loc='center')
        ax1.axis('off')

        fig.savefig(
            f"~/ch5_integration/5.4 graph_localisation/5.4.1 clusters/plots/combo_gwas_proteomics_enrichment_cooccurrence_{analysis_category_with_genegwas}.png",
            bbox_inches="tight", dpi=100)
        fig.savefig(
            f"~/ch5_integration/5.4 graph_localisation/5.4.1 clusters/plots/combo_gwas_proteomics_enrichment_cooccurrence_{analysis_category_with_genegwas}.svg",
            bbox_inches="tight", dpi=100)

        plt.close()
