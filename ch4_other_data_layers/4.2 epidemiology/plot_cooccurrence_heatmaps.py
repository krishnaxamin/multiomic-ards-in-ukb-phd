""" Plot alpha-MLE for the 3 cohorts (genomics, proteomics, prior proteomics) as heatmaps. """
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from utils.significance_labelling import significance_labelling

mpl.use('Agg')

""" Plotting function """


def plot_heatmap(input_data: pd.DataFrame, omic: str):
    data_to_plot = input_data[['disease1', 'disease2', 'alpha_hat']].copy()

    data_to_plot = (data_to_plot
                    .merge(disease_info[['disease_field', 'code_chapter']]
                           .rename(columns={'disease_field': 'disease1'}))
                    .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                    .merge(disease_info[['disease_field', 'code_chapter']]
                           .rename(columns={'disease_field': 'disease2'}))
                    .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))

    data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                              categories=sorted(data_to_plot['disease1'].unique()))
    data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                              categories=sorted(data_to_plot['disease2'].unique()))

    # map data to symmetric categorical index space
    all_labels = sorted(set(data_to_plot['disease1']) | set(data_to_plot['disease2']))
    # make mapping have 1-tick gaps between chapters, so that there are whitespace gaps/boundaries between chapters
    label_to_num = {}
    current_pos = 0
    for i in range(len(all_labels)):
        if i == 0:
            label_to_num[all_labels[i]] = current_pos
            current_pos += 1
            continue
        current_chapter = int(re.search(r"\((\d+)\)", all_labels[i]).group(1))
        past_chapter = int(re.search(r"\((\d+)\)", all_labels[i - 1]).group(1))
        if past_chapter < current_chapter:
            current_pos += 1  # insert a 1-tick wide gap into the axis at the boundary between chapters
        label_to_num[all_labels[i]] = current_pos
        current_pos += 1
    # label_to_num = {lab: i for i, lab in enumerate(all_labels)}
    data_to_plot['disease1_as_num'] = data_to_plot['disease1'].map(label_to_num)
    data_to_plot['disease2_as_num'] = data_to_plot['disease2'].map(label_to_num)

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.scatterplot(data_to_plot, x='disease1_as_num', y='disease2_as_num', hue='alpha_hat',
                    palette='viridis', zorder=3, ax=ax, s=50)
    # plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', title='Co-occurrence (alpha)', fontsize=8)
    plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=9)
    plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=9)
    plt.grid(True, zorder=0)

    if ax.yaxis_inverted():
        ax.invert_yaxis()

    # set ticks labels to the string versions
    tick_positions = [label_to_num[label] for label in all_labels]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(all_labels, rotation=90, size=8)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(all_labels, size=8)

    # add diagonal line going top-left to bottom-right, ensuring whitespaces are dealt with properly
    coords = np.array(tick_positions)
    ax.plot(coords, coords,
            color='black', linestyle='--', linewidth=0.5)

    # build colorbar
    norm = mpl.colors.Normalize(vmin=data_to_plot.alpha_hat.min(),
                                vmax=data_to_plot.alpha_hat.max())
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])  # required boilerplate

    # remove legend, attach colorbar
    ax.get_legend().remove()
    plt.colorbar(sm, ax=ax, label='Co-occurrence (alpha)')

    plt.tight_layout()
    plt.savefig(f"~/ch4_other_data_layers/4.2 epidemiology/plots/heatmap_alphamle_{omic}.svg")
    plt.savefig(f"~/ch4_other_data_layers/4.2 epidemiology/plots/heatmap_alphamle_{omic}.png")
    plt.close()


""" Read in data """
genomics = pd.read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
proteomics = pd.read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
prior_proteomics = pd.read_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences_alphamle.csv')

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
code_chapter_to_field_mapper = dict(zip(disease_info.code_chapter, disease_info.disease_field))

""" FDR data """
genomics = significance_labelling(genomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                  pair_columns=['disease1', 'disease2'])
genomics_sig = genomics[genomics['fdr_sig'] == 1].copy()
proteomics = significance_labelling(proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                    pair_columns=['disease1', 'disease2'])
proteomics_sig = proteomics[proteomics['fdr_sig'] == 1].copy()
prior_proteomics = significance_labelling(prior_proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                          pair_columns=['disease1', 'disease2'])
prior_proteomics_sig = prior_proteomics[prior_proteomics['fdr_sig'] == 1].copy()

""" Plot """
plot_heatmap(genomics_sig, 'genomics')
plot_heatmap(proteomics_sig, 'proteomics')
plot_heatmap(prior_proteomics_sig, 'prior_proteomics')
