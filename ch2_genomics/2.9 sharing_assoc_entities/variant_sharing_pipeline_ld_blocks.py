""" Pipeline to perform between-disease variant sharing analyses using LD blocks: Fisher tests, cosine similarities.
This time, commonalities calculated between LD block hits, where a pre-computed LD block (from Berisa and Pickrell 2016)
is a hit if it has an assoc variant within it (method from Donertas et al, 2021). """
from pandas import read_csv, DataFrame, concat, Categorical, Series
from collections import defaultdict
from scipy.stats import spearmanr, fisher_exact
from itertools import combinations
from utils.significance_labelling import significance_labelling

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import sklearn.metrics
import pyprind
import sys

mpl.use('TkAgg')

""" General info """
disease_info = read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
ld_blocks = read_csv('~/data/external/genomics/berisa_pickrell_ld_blocks.bed', sep='\s+')
ld_blocks['block'] = ['block' + str(x) for x in ld_blocks.index]
ld_blocks['chr'] = ld_blocks['chr'].apply(lambda x: x.replace('chr', ''))

""" Collate variants """
info_dict = {'variant_type': 'all-types',
             'assoc_label': 'assoc',
             'info_label': 'high-info'}
file_suffix = '_'.join([
    info_dict['variant_type'],
    info_dict['assoc_label'],
    info_dict['info_label']])

variants = read_csv('ukbiobank/gwas/saige/pan-ukbb-eur_assoc_saige_lifestyles.csv')
variants['uniqueID'] = variants['CHR'].astype(str) + '_' + variants['POS'].astype(str) + '_' + variants[
    'UKB_ref'] + '_' + variants['UKB_alt']
variants = variants.drop_duplicates(['uniqueID', 'disease'], ignore_index=True).astype({'CHR': 'str'})
variants.loc[variants['Allele2_is_UKB_alt'] == 0, 'Allele2_is_UKB_alt'] = -1

# filter variants
variants = variants[(variants['assoc_label'] == 'assoc') & (variants['info_label'] == 'high_info')].copy()

""" Convert variant hits to LD block hits """


def match_variant_to_ld_block(variant_row, ld_blocks_df=ld_blocks):
    try:
        matched_block = ld_blocks_df[(ld_blocks_df['chr'] == variant_row['CHR']) &
                                     (ld_blocks_df['start'] < int(variant_row['POS'])) &
                                     (ld_blocks_df['stop'] > int(variant_row['POS']))]['block'].values[0]
    except IndexError:
        matched_block = 'no_match'

    return matched_block


def variant_hits_to_ld_block_hits(variant_set_to_convert):
    variants_with_ld_blocks = variant_set_to_convert.copy()
    variants_with_ld_blocks['ld_block'] = variants_with_ld_blocks.apply(match_variant_to_ld_block, axis=1)
    variants_with_ld_blocks = variants_with_ld_blocks[variants_with_ld_blocks['ld_block'] != 'no_match'].copy()

    ld_block_hits = DataFrame()
    for disease in variants_with_ld_blocks.disease.unique():
        disease_variants = variants_with_ld_blocks[variants_with_ld_blocks['disease'] == disease][
            ['BETA', 'p.value', 'Allele2_is_UKB_alt', 'uniqueID', 'disease', 'ld_block']].drop_duplicates(
            ['uniqueID', 'disease', 'ld_block'], ignore_index=True).copy()
        for block in disease_variants.ld_block.unique():
            disease_variants_block = disease_variants[disease_variants['ld_block'] == block].copy()
            block_index_beta = disease_variants_block.sort_values(by='p.value', ascending=True)['BETA'].values[0]
            disease_variants_block.loc[disease_variants_block['p.value'] == 0, 'p.value'] = \
                disease_variants_block[disease_variants_block['p.value'] != 0]['p.value'].min()
            block_weighted_avg_beta = (disease_variants_block['BETA'] * (
                    -np.log10(disease_variants_block['p.value']) / -np.log10(
                disease_variants_block['p.value']).sum())).sum() / len(disease_variants_block)
            ld_block_hits = concat([ld_block_hits, DataFrame([{'disease': disease,
                                                               'ld_block': block,
                                                               'index_beta': block_index_beta,
                                                               'weighted_avg_beta': block_weighted_avg_beta}])])

    return ld_block_hits.reset_index(drop=True)


ld_block_hits = variant_hits_to_ld_block_hits(variant_set_to_convert=variants)

""" Get number of LD blocks shared """


def get_ld_blocks_shared(ld_block_set_to_explore, file_suffix):
    # digest output_file_info_dict
    diseases_with_shared_ld_blocks_file_name = 'ld_block_sharing_' + file_suffix

    # identify possible disease pairs
    disease_pairs = list(combinations(sorted(list(ld_block_set_to_explore.disease.unique())), 2))

    # count number of shared LD blocks per disease pair
    diseases_with_shared_ld_blocks_df = DataFrame()
    for disease_pair in disease_pairs:
        shared_ld_blocks = (
                set(ld_block_set_to_explore[ld_block_set_to_explore['disease'] == disease_pair[0]]['ld_block']) &
                set(ld_block_set_to_explore[ld_block_set_to_explore['disease'] == disease_pair[1]]['ld_block']))
        if len(shared_ld_blocks) > 0:
            diseases_with_shared_ld_blocks_df = concat(
                [diseases_with_shared_ld_blocks_df, DataFrame([{'disease1': disease_pair[0],
                                                                'disease2': disease_pair[1],
                                                                'count': len(shared_ld_blocks)}])])

    return diseases_with_shared_ld_blocks_df.reset_index(drop=True)


diseases_with_shared_blocks_df = get_ld_blocks_shared(ld_block_hits, file_suffix=file_suffix)

""" Find (significant) LD block sharing """
# Firth co-occurrence tests
num_blocks = len(ld_blocks)


# function for a single row-as-dict
def fisher_testing(shared_blocks_data_row):
    disease1 = shared_blocks_data_row['disease1']
    disease2 = shared_blocks_data_row['disease2']

    num_blocks_disease1 = len(ld_block_hits[ld_block_hits['disease'] == disease1]['ld_block'].unique())
    num_blocks_disease2 = len(ld_block_hits[ld_block_hits['disease'] == disease2]['ld_block'].unique())

    num_shared_blocks = shared_blocks_data_row['count']

    a = num_shared_blocks
    b = num_blocks_disease2 - num_shared_blocks
    c = num_blocks_disease1 - num_shared_blocks
    d = num_blocks - num_blocks_disease1 - b

    fisher_result = fisher_exact([[a, b], [c, d]], alternative='greater')

    return {'disease1': disease1, 'disease2': disease2, 'n_shared': num_shared_blocks,
            'odds_ratio': fisher_result.statistic, 'pval': fisher_result.pvalue}


fisher_dicts = []
diseases_with_shared_blocks_df_as_dicts = diseases_with_shared_blocks_df.to_dict(orient='records')
for row_dict in diseases_with_shared_blocks_df_as_dicts:
    fisher_dicts.append(fisher_testing(row_dict))
fisher_tests_df = DataFrame.from_records(fisher_dicts)
fisher_tests_df_sig_labelled = significance_labelling(fisher_tests_df)

""" Calculate effect size similarity between LD blocks shared between disease pairs """
variants_ld_blocks = variants.copy()
variants_ld_blocks['ld_block'] = variants_ld_blocks.apply(match_variant_to_ld_block, axis=1)
variants_ld_blocks = variants_ld_blocks[variants_ld_blocks['ld_block'] != 'no_match'].copy()

disease_pair_list = [(diseases_with_shared_blocks_df['disease1'][i], diseases_with_shared_blocks_df['disease2'][i]) for
                     i in range(len(diseases_with_shared_blocks_df))]

# results store
ld_block_correlations = defaultdict()

bar = pyprind.ProgBar(len(disease_pair_list), stream=sys.stdout, title='Calculating cossims')
for disease_pair in disease_pair_list:
    # disease_pair = ['J12', 'M81']
    ld_block_correlations[tuple(disease_pair)] = defaultdict(dict)
    shared_ld_blocks = (set(variants_ld_blocks[variants_ld_blocks['disease'] == disease_pair[0]]['ld_block']) &
                        set(variants_ld_blocks[variants_ld_blocks['disease'] == disease_pair[1]]['ld_block']))
    variants_shared_ld_blocks = variants_ld_blocks[(variants_ld_blocks['disease'].isin(disease_pair)) &
                                                   (variants_ld_blocks['ld_block'].isin(list(shared_ld_blocks)))].copy()

    # take the range of POS in variants_shared_ld_blocks, split into 10kb boxes and for each disease, collate the mean beta for that box
    variants_shared_ld_blocks_dfs = [d for block, d in variants_shared_ld_blocks.groupby('ld_block')]
    for block, df in variants_shared_ld_blocks.groupby('ld_block'):
        bin_size = 10000
        pos_range = [df.POS.min(), df.POS.max()]
        if df.POS.min() == df.POS.max():
            pos_range[1] += 1
        pos_bins = [(i, min(i + bin_size, pos_range[1])) for i in range(pos_range[0], pos_range[1], bin_size)]
        disease1_betas = []
        disease2_betas = []
        for pos_bin in pos_bins:
            disease1_betas.append(
                df[(df['POS'] >= pos_bin[0]) & (df['POS'] <= pos_bin[1]) & (df['disease'] == disease_pair[0])][
                    'BETA'].mean())
            disease2_betas.append(
                df[(df['POS'] >= pos_bin[0]) & (df['POS'] <= pos_bin[1]) & (df['disease'] == disease_pair[1])][
                    'BETA'].mean())
        disease1_betas = [0 if np.isnan(x) else x for x in disease1_betas]
        disease2_betas = [0 if np.isnan(x) else x for x in disease2_betas]
        # correlate the mean betas across the 10kb boxes
        ld_block_cossim = sklearn.metrics.pairwise.cosine_similarity([disease1_betas, disease2_betas])
        ld_block_correlations[tuple(disease_pair)][block] = ld_block_cossim[0, 1]
    bar.update()

# export as df
ld_block_correlations_df = DataFrame()
for (d1, d2), ld_block_dict in ld_block_correlations.items():
    for block, cossim in ld_block_dict.items():
        ld_block_correlations_df = concat([ld_block_correlations_df, DataFrame([{'disease1': d1,
                                                                                 'disease2': d2,
                                                                                 'ld_block': block,
                                                                                 'cossim': cossim}])])

ld_block_correlations_df.to_csv(
    "~/data/internal/genomics/shared_assoc_entities/ld_block_blockwise_cossims.csv",
    index=False)

""" Plot combined heatmap for LD block overlap and correlations """
# colour = average corr, size = num LD blocks, shape = FDR sig status

# apply groupby and aggregation
ld_blocks_cossim_agg = ld_block_correlations_df.groupby(['disease1', 'disease2'])['cossim'].mean().reset_index()
# combine
ld_blocks_combined_info = ld_blocks_cossim_agg.merge(
    fisher_tests_df_sig_labelled[['disease1', 'disease2', 'n_shared', 'fdr_sig']])
ld_blocks_combined_info.to_csv(
    "~/data/internal/genomics/shared_assoc_entities/ld_block_sharing_fisher_cossims.csv",
    index=False)
ld_blocks_for_heatmap = concat(
    [ld_blocks_combined_info, ld_blocks_combined_info.rename(columns={'disease1': 'disease2', 'disease2': 'disease1'})])

# prep data for plotting
ld_blocks_for_heatmap = (ld_blocks_for_heatmap
                         .merge(disease_info[['icd10_three_letter', 'code_chapter']]
                                .rename(columns={'icd10_three_letter': 'disease1'}))
                         .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                         .merge(disease_info[['icd10_three_letter', 'code_chapter']]
                                .rename(columns={'icd10_three_letter': 'disease2'}))
                         .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
ld_blocks_for_heatmap['disease1'] = Categorical(ld_blocks_for_heatmap['disease1'],
                                                categories=sorted(ld_blocks_for_heatmap['disease1'].unique()))
ld_blocks_for_heatmap['disease2'] = Categorical(ld_blocks_for_heatmap['disease2'],
                                                categories=sorted(ld_blocks_for_heatmap['disease2'].unique()))
ld_blocks_for_heatmap['fdr_sig'] = ld_blocks_for_heatmap['fdr_sig'].astype(int).astype(str)
ld_blocks_for_heatmap.loc[ld_blocks_for_heatmap['fdr_sig'] == '1', 'fdr_sig'] = 'Yes'
ld_blocks_for_heatmap.loc[ld_blocks_for_heatmap['fdr_sig'] == '0', 'fdr_sig'] = 'No'

norm = mpl.colors.TwoSlopeNorm(
    vcenter=0,
    vmin=ld_blocks_for_heatmap['cossim'].min(),
    vmax=ld_blocks_for_heatmap['cossim'].max()
)

# map data to symmetric categorical index space
all_labels = sorted(set(ld_blocks_for_heatmap['disease1']) | set(ld_blocks_for_heatmap['disease2']))
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
ld_blocks_for_heatmap['disease1_as_num'] = ld_blocks_for_heatmap['disease1'].map(label_to_num)
ld_blocks_for_heatmap['disease2_as_num'] = ld_blocks_for_heatmap['disease2'].map(label_to_num)

# plot
fig, ax = plt.subplots(figsize=(12, 10))
sns.scatterplot(ld_blocks_for_heatmap.rename(columns={'fdr_sig': 'FDR < 0.05',
                                                      'n_shared': 'Number LD \nblocks shared',
                                                      'cossim': 'Mean cosine \nsimilarity'}),
                x='disease1_as_num', y='disease2_as_num', hue='Mean cosine \nsimilarity',
                size='Number LD \nblocks shared', style='FDR < 0.05', palette=plt.cm.coolwarm, hue_norm=norm,
                style_order=['Yes', 'No'], zorder=3, edgecolor=None, sizes=(50, 150),
                markers={'Yes': 'o', 'No': '^'}, ax=ax)
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('Disease ICD-10 code (ICD-10 chapter)')
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
    if handle._label in ['FDR < 0.05', 'Number LD \nblocks shared', 'Mean cosine \nsimilarity']:  # section header
        text.set_fontweight('bold')

plt.tight_layout()
plt.savefig(f"~/ch2_genomics/2.9 sharing_assoc_entities/plots/ld_block_sharing_fisher_cossims.png")
plt.savefig(f"~/ch2_genomics/2.9 sharing_assoc_entities/plots/ld_block_sharing_fisher_cossims.svg")
plt.close()
