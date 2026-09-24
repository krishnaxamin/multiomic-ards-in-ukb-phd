""" Script to analyse sharing of associated proteins, plot results, and conduct auxiliary analysis and plotting. """
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
import collections
import sys
import scipy
import sklearn.metrics
import re
import numpy as np

from pyprind import ProgBar
from utils.significance_labelling import significance_labelling

mpl.use('TkAgg')

""" Read in data """
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'


def analyse_sharing(dataset1: pd.DataFrame, dataset2: pd.DataFrame, export_path: str):
    """
    Test sharing significance and cosine similarity between effect sizes, given two datasets of protein-disease
    associations.
    :param dataset1:
    :param dataset2:
    :param export_path:
    :return:
    """

    # dataset1
    dataset1_proteins_per_disease = {disease: set(df.protein) for disease, df in dataset1.groupby('disease')}

    # dataset2
    dataset2_proteins_per_disease = {disease: set(df.protein) for disease, df in dataset2.groupby('disease')}

    disease_pairs = list(itertools.combinations(sorted(list(set(dataset1.disease) & set(dataset2.disease))), 2))

    shared_protein_dict = collections.defaultdict(list)
    sharing_results = []
    bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Identifying shared proteins')
    for disease_pair in disease_pairs:
        disease1 = disease_pair[0]
        disease2 = disease_pair[1]

        # get shared proteins
        shared_protein = list(dataset1_proteins_per_disease[disease1] & dataset2_proteins_per_disease[disease2])
        shared_protein_dict[disease_pair] = shared_protein

        # test sharing significance
        num_shared = len(shared_protein)
        num_assoc_disease1 = len(dataset1_proteins_per_disease[disease1])
        num_assoc_disease2 = len(dataset2_proteins_per_disease[disease2])

        # set up contingency table [[a, b], [c, d]] where columns are in-test-set +/-, rows are annotation +/-
        #                    disease 1
        #                -------------------
        #                |   +   |   -     |
        #  --------------|-------|---------|
        #   disease 2 +  |   a   |   b     |  num_assoc_disease2
        #   disease 2 -  |   c   |   d     |
        #  --------------|-------|---------|
        #                 num_assoc_disease1  n_tested_both

        a = num_shared
        b = num_assoc_disease2 - num_shared
        c = num_assoc_disease1 - num_shared
        d = 2920 - num_assoc_disease1 - b

        fisher_result = scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        # cosine similarity
        if num_shared > 0:
            shared_protein_effects_paired = pd.merge(
                dataset1[dataset1.disease == disease1].drop(columns='disease').rename(
                    columns={'beta': 'disease1_beta'}),
                dataset2[dataset2.disease == disease2].drop(columns='disease').rename(
                    columns={'beta': 'disease2_beta'}),
            )
            cossim = \
                sklearn.metrics.pairwise.cosine_similarity(
                    shared_protein_effects_paired[['disease1_beta', 'disease2_beta']].T)[
                    0, 1]
        else:
            cossim = 0

        sharing_results.append({'disease1': disease1, 'disease2': disease2, 'n_shared': num_shared,
                                'odds_ratio': fisher_result.statistic, 'pval': fisher_result.pvalue,
                                'cossim': cossim})

        bar.update()

    sharing_results_df = pd.DataFrame(sharing_results)
    sharing_results_df = sharing_results_df[sharing_results_df.n_shared > 0].copy()
    sharing_results_df = significance_labelling(sharing_results_df)
    sharing_results_df.to_csv(f"{export_path}.csv", index=False)

    return sharing_results_df


""" Plot """


# as in variant and gene sharing
# circle = FDR < 0.05; triangle = FDR insigh
# colour = cossim
# size = n_proteins shared


def plot_asymmetric_heatmap(sharing_df1: pd.DataFrame, sharing_df2: pd.DataFrame, plot_path: str,
                            plot_xaxis_label: str = 'Disease ICD-10 code (ICD-10 chapter)',
                            plot_yaxis_label: str = 'Disease ICD-10 code (ICD-10 chapter)'):
    """
    Given two different datasets describing protein sharing, plot both on the same set of axes in a scatterplot-based
    heatmap. style=FDR-sig; hue=cossim; size=log10(n_proteins_shared)
    :param sharing_df1:
    :param sharing_df2:
    :param plot_path:
    :param plot_xaxis_label:
    :param plot_yaxis_label:
    :return:
    """

    data_to_plot = pd.concat([sharing_df1[['disease1', 'disease2', 'n_shared', 'fdr_sig', 'cossim']],
                              sharing_df2[['disease1', 'disease2', 'n_shared', 'fdr_sig', 'cossim']]
                             .rename(columns={'disease1': 'disease2', 'disease2': 'disease1'})])

    data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                              categories=sorted(data_to_plot['disease1'].unique()))
    data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                              categories=sorted(data_to_plot['disease2'].unique()))
    data_to_plot['fdr_sig'] = data_to_plot['fdr_sig'].astype(int).astype(str)
    data_to_plot.loc[data_to_plot['fdr_sig'] == '1', 'fdr_sig'] = 'Yes'
    data_to_plot.loc[data_to_plot['fdr_sig'] == '0', 'fdr_sig'] = 'No'

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

    norm = mpl.colors.TwoSlopeNorm(
        vcenter=0,
        vmin=data_to_plot['cossim'].min(),
        vmax=data_to_plot['cossim'].max()
    )

    # plot
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.scatterplot(data=data_to_plot.assign(log_n_shared=np.log10(data_to_plot['n_shared'])).rename(
        columns={'fdr_sig': 'FDR < 0.05', 'log_n_shared': 'log10(number \nproteins \nshared)',
                 'cossim': 'Cosine \nsimilarity'}
    ),
        x='disease1_as_num', y='disease2_as_num', hue='Cosine \nsimilarity', style='FDR < 0.05',
        style_order=['Yes', 'No'],
        markers={'Yes': 'o', 'No': '^'}, edgecolor=None, palette=plt.cm.coolwarm, hue_norm=norm,
        size='log10(number \nproteins \nshared)', sizes=(25, 160),
        zorder=3, ax=ax)
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.xlabel(plot_xaxis_label, size=11)
    plt.ylabel(plot_yaxis_label, size=11)
    plt.grid(True, zorder=0)

    if ax.yaxis_inverted():
        ax.invert_yaxis()

    # set ticks labels to the string versions
    tick_positions = [label_to_num[label] for label in all_labels]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(all_labels, rotation=90, size=9)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(all_labels, size=9)

    # add diagonal line going top-left to bottom-right, ensuring whitespaces are dealt with properly
    coords = np.array(tick_positions)
    ax.plot(coords, coords,
            color='black', linestyle='--', linewidth=0.5)

    # bold legend section titles
    legend = ax.get_legend()
    for text, handle in zip(legend.texts, legend.legend_handles):
        if handle._label in ['FDR < 0.05', 'log10(number \nproteins \nshared)',
                             'Cosine \nsimilarity']:  # section header
            text.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(f"{plot_path}.png")
    plt.savefig(f"{plot_path}.svg")
    plt.close()

    return


""" Firth-Firth and Cox-Cox """
# firth-firth = uppper; cox-cox = lower
firth = pd.read_csv("~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_firth.csv")
sig_firth = firth[firth['fdr_sig'] == 1][['disease', 'protein', 'beta']].copy()
sig_firth = sig_firth.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(
    columns='disease').rename(columns={'code_chapter': 'disease'})

cox = pd.read_csv("~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_cox.csv")
sig_cox = cox[cox['fdr_sig'] == 1][['disease', 'protein', 'beta']].copy()
sig_cox = sig_cox.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(
    columns='disease').rename(columns={'code_chapter': 'disease'})

firth_sharing = analyse_sharing(dataset1=sig_firth, dataset2=sig_firth,
                                export_path=f"~/data/internal/proteomics/firth_protein_sharing_fisher_cossims")
cox_sharing = analyse_sharing(dataset1=sig_cox, dataset2=sig_cox,
                              export_path=f"~/data/internal/proteomics/cox_protein_sharing_fisher_cossims")

plot_asymmetric_heatmap(sharing_df1=firth_sharing, sharing_df2=cox_sharing,
                        plot_path=f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/firthfirth_coxcox_protein_sharing_fisher_cossims")

""" Firth-Cox and Cox-Firth """
firth = pd.read_csv(
    f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_firth.csv")
sig_firth = firth[firth['fdr_sig'] == 1][['disease', 'protein', 'beta']].copy()
sig_firth = sig_firth.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(
    columns='disease').rename(columns={'code_chapter': 'disease'})

cox = pd.read_csv(
    f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_cox.csv")
sig_cox = cox[cox['fdr_sig'] == 1][['disease', 'protein', 'beta']].copy()
sig_cox = sig_cox.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})).drop(
    columns='disease').rename(columns={'code_chapter': 'disease'})

firth_cox_sharing = analyse_sharing(dataset1=sig_firth, dataset2=sig_cox,
                                    export_path=f"~/data/internal/proteomics/shared_assoc_proteins/firth_cox_protein_sharing_fisher_cossims")
cox_firth_sharing = analyse_sharing(dataset1=sig_cox, dataset2=sig_firth,
                                    export_path=f"~/data/internal/proteomics/shared_assoc_proteins/cox_firth_protein_sharing_fisher_cossims")

plot_asymmetric_heatmap(sharing_df1=firth_cox_sharing, sharing_df2=cox_firth_sharing,
                        plot_path=f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/firthcox_coxfirth_protein_sharing_fisher_cossims",
                        plot_xaxis_label='Disease ICD-10 code (ICD-10 chapter) - Firth',
                        plot_yaxis_label='Disease ICD-10 code (ICD-10 chapter) - Cox')

""" Auxiliary analyses """
firth_sharing_chapter_overlap = firth_sharing[firth_sharing.fdr_sig == 1].merge(
    disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease1', 'icd10_chapter': 'disease1_chapter'})
).merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease2', 'icd10_chapter': 'disease2_chapter'}))
firth_sharing_chapter_overlap['chapter_match'] = firth_sharing_chapter_overlap['disease1_chapter'] == firth_sharing_chapter_overlap['disease2_chapter']
firth_sharing_chapter_overlap['chapter_match'].sum()
# enrichment of within-chapter pairs
X = firth_sharing_chapter_overlap['chapter_match'].sum()
Y = len(firth_sharing_chapter_overlap)
scipy.stats.fisher_exact([[X, 293 - X], [Y - X, 1985 - Y + X]], alternative='greater')

cox_sharing_chapter_overlap = cox_sharing[cox_sharing.fdr_sig == 1].merge(
    disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease1', 'icd10_chapter': 'disease1_chapter'})
).merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease2', 'icd10_chapter': 'disease2_chapter'}))
cox_sharing_chapter_overlap['chapter_match'] = cox_sharing_chapter_overlap['disease1_chapter'] == cox_sharing_chapter_overlap['disease2_chapter']
cox_sharing_chapter_overlap['chapter_match'].sum()
# enrichment of within-chapter pairs
X = cox_sharing_chapter_overlap['chapter_match'].sum()
Y = len(cox_sharing_chapter_overlap)
scipy.stats.fisher_exact([[X, 293 - X], [Y - X, 1985 - Y + X]], alternative='greater')

firth_cox_sharing_chapter_overlap = pd.concat([firth_cox_sharing[firth_cox_sharing.fdr_sig == 1],
                                               cox_firth_sharing[cox_firth_sharing.fdr_sig == 1]]).merge(
    disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease1', 'icd10_chapter': 'disease1_chapter'})
).merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease2', 'icd10_chapter': 'disease2_chapter'}))
firth_cox_sharing_chapter_overlap['chapter_match'] = firth_cox_sharing_chapter_overlap['disease1_chapter'] == firth_cox_sharing_chapter_overlap['disease2_chapter']
firth_cox_sharing_chapter_overlap['chapter_match'].sum()
# enrichment of within-chapter pairs
X = firth_cox_sharing_chapter_overlap['chapter_match'].sum()
Y = len(firth_cox_sharing_chapter_overlap)
scipy.stats.fisher_exact([[X, 293*2 - X], [Y - X, 1985*2 - Y + X]], alternative='greater')

""" Auxiliary plots """
# diseases participating most in significant protein sharing from Firth regressions
firth_sharing_freq = pd.Series(firth_sharing[firth_sharing.fdr_sig == 1].disease1.to_list() + firth_sharing[firth_sharing.fdr_sig == 1].disease2.to_list()).value_counts()
plt.figure(figsize=(8, 6))
sns.barplot(pd.DataFrame(firth_sharing_freq).reset_index().set_axis(['disease', 'count'], axis=1).sort_values(by='count', ascending=False),
            x='disease', y='count')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('Number of FDR < 0.05 ARD pairs involving this ARD')
plt.gca().tick_params(axis='x', rotation=90)
for container in plt.gca().containers:
    plt.gca().bar_label(container)
plt.grid(True, which='major', axis='y', color='lightgrey', lw=0.8, ls='--')
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig(f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/firth_protein_sharing_fisher_disease_freq.png")
plt.savefig(f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/firth_protein_sharing_fisher_disease_freq.svg")
plt.close()

# diseases participating most in dissimilar significant protein sharing from Cox regressions
neg_cox_sharing_freq = pd.Series(list(cox_sharing[(cox_sharing.fdr_sig == 1) & (cox_sharing.cossim < 0)].disease1) + list(cox_sharing[(cox_sharing.fdr_sig == 1) & (cox_sharing.cossim < 0)].disease2)).value_counts()
plt.figure(figsize=(8, 6))
sns.barplot(pd.DataFrame(neg_cox_sharing_freq).reset_index().set_axis(['disease', 'count'], axis=1).sort_values(by='count', ascending=False),
            x='disease', y='count')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('Number of dissimilar FDR < 0.05 ARD pairs involving this ARD')
plt.gca().tick_params(axis='x', rotation=90)
for container in plt.gca().containers:
    plt.gca().bar_label(container)
plt.grid(True, which='major', axis='y', color='lightgrey', lw=0.8, ls='--')
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig(f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/cox_protein_sharing_fisher_dissimilar_disease_freq.png")
plt.savefig(f"~/ch3_proteomics/3.5 sharing_assoc_proteins/plots/cox_protein_sharing_fisher_dissimilar_disease_freq.svg")
plt.close()
