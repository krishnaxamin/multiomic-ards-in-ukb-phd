"""
Identify high-confidence colocalised variants and map these to Ensembl hg19 genes.
Plot the colocalisation with gene mapping.
"""
from typing import Dict

import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('TkAgg')


# for both abf and susie, identify the 4 categories of colocalisation
# 1. strong coloc (PP.H4 > 0.9), strong variant (SNP.PP.H4 > 0.9)
# 2. strong coloc (PP.H4 > 0.9), moderate variant (0.9 >= SNP.PP.H4 > 0.5)
# 3. moderate coloc (0.9 >= PP.H4 > 0.5), strong variant (SNP.PP.H4 > 0.9)
# 4. moderate coloc (0.9 >= PP.H4 > 0.5), moderate variant (0.9 >= SNP.PP.H4 > 0.5)


def abf_find_colocs(abf_summary: pd.DataFrame, strong_threshold: float = 0.9, moderate_threshold: float = 0.5) -> Dict[
    str, pd.DataFrame]:
    """
    Find regions where colocalisation is likely, from single-causal-variant ABF results.
    :param abf_summary:
    :param strong_threshold:
    :param moderate_threshold:
    :return:
    """
    strong_pairs = abf_summary[abf_summary['PP.H4.abf'] > strong_threshold].copy()
    moderate_pairs = abf_summary[(abf_summary['PP.H4.abf'] <= strong_threshold) &
                                 (abf_summary['PP.H4.abf'] > moderate_threshold)].copy()

    return {'strong_coloc': strong_pairs[['combo', 'region', 'PP.H4.abf']],
            'moderate_coloc': moderate_pairs[['combo', 'region', 'PP.H4.abf']]}


def abf_probable_variant(
        abf_results: pd.DataFrame, coloc_dict: Dict[str, pd.DataFrame], strong_threshold: float = 0.9,
        moderate_threshold: float = 0.5) -> pd.DataFrame:
    """
    Identify the likely causal variant behind colocalisations, from single-causal-variant ABF results.
    :param abf_results:
    :param coloc_dict:
    :param strong_threshold:
    :param moderate_threshold:
    :return:
    """
    abf_results_working = abf_results[['snp', 'SNP.PP.H4', 'combo', 'region']].copy()
    results_df_list = []
    for cat, cat_df in coloc_dict.items():
        for combo_region, combo_region_df in cat_df.groupby(['combo', 'region']):
            assert len(combo_region_df) == 1
            variant_df = abf_results_working.merge(combo_region_df)
            strong_variants = variant_df[variant_df['SNP.PP.H4'] > strong_threshold].assign(
                category=f"{cat}-strong_variant")
            moderate_variants = variant_df[(variant_df['SNP.PP.H4'] <= strong_threshold) &
                                           (variant_df['SNP.PP.H4'] > moderate_threshold)].assign(
                category=f"{cat}-moderate_variant")
            results_df_list.append(pd.concat([strong_variants, moderate_variants]))

    # snp, combo, region, PP.H4.abf, SNP.PP.H4, category (one of the 4 defined above)
    return pd.concat(results_df_list)


def susie_find_colocs(susie_summary: pd.DataFrame, strong_threshold: float = 0.9, moderate_threshold: float = 0.5) -> \
        Dict[str, pd.DataFrame]:
    """
    Find regions where colocalisation is likely, from multiple-causal-variant SuSiE results.
    :param susie_summary:
    :param strong_threshold:
    :param moderate_threshold:
    :return:
    """
    susie_summary_working = susie_summary[['combo', 'region', 'PP.H4.abf']].copy()

    # add a counter of how many rows there are for each combo-region, for reference with the 'results' df
    g = susie_summary.groupby(['combo', 'region'])
    sizes = g['combo'].transform('size')
    counts = g.cumcount()
    susie_summary_working['pp_h4_suffix'] = np.where(
        sizes == 1,
        'abf',
        'row' + (counts + 1).astype(str)
    )

    strong_pairs = susie_summary_working[susie_summary_working['PP.H4.abf'] > strong_threshold].copy()
    moderate_pairs = susie_summary_working[(susie_summary_working['PP.H4.abf'] <= strong_threshold) &
                                           (susie_summary_working['PP.H4.abf'] > moderate_threshold)].copy()

    return {'strong_coloc': strong_pairs, 'moderate_coloc': moderate_pairs}


def susie_probable_variant(
        susie_results: pd.DataFrame, coloc_dict: Dict[str, pd.DataFrame], strong_threshold: float = 0.9,
        moderate_threshold: float = 0.5) -> pd.DataFrame:
    """
    Identify the likely causal variant behind colocalisations, from multiple-causal-variant SuSiE results.
    :param susie_results:
    :param coloc_dict:
    :param strong_threshold:
    :param moderate_threshold:
    :return:
    """
    susie_results_working = susie_results.copy()
    results_df_list = []
    for cat, cat_df in coloc_dict.items():
        for combo_region, combo_region_df in cat_df.groupby(['combo', 'region']):
            # assert len(combo_region_df) == 1, f"{combo_region_df}"
            for _, row in combo_region_df.iterrows():
                pp_h4_suffix = row['pp_h4_suffix']
                variant_df = susie_results_working.merge(pd.DataFrame([row]))
                strong_variants = variant_df[variant_df[f"SNP.PP.H4.{pp_h4_suffix}"] > strong_threshold].assign(
                    category=f"{cat}-strong_variant")
                moderate_variants = variant_df[(variant_df[f"SNP.PP.H4.{pp_h4_suffix}"] <= strong_threshold) &
                                               (variant_df[f"SNP.PP.H4.{pp_h4_suffix}"] > moderate_threshold)].assign(
                    category=f"{cat}-moderate_variant")
                results_df_list.append(pd.concat([strong_variants, moderate_variants]).rename(
                    columns={f"SNP.PP.H4.{pp_h4_suffix}": 'SNP.PP.H4'})[
                                           ['snp', 'PP.H4.abf', 'SNP.PP.H4', 'combo', 'region', 'category']])

    # snp, combo, region, PP.H4.susie, SNP.PP.H4, category (one of the 4 defined above)
    return pd.concat(results_df_list)


abf_summaries = pd.read_csv(
    f"~/ch2_genomics/2.13 colocalisations/results/coloc_summaries_abf.csv")
abf_results = pd.read_csv(f"~/ch2_genomics/2.13 colocalisations/results/coloc_results_abf.csv")

abf_processed = abf_probable_variant(abf_results=abf_results,
                                     coloc_dict=abf_find_colocs(abf_summary=abf_summaries))

susie_summaries = pd.read_csv(
    f"~/ch2_genomics/2.13 colocalisations/results/coloc_summaries_susie.csv")
susie_results = pd.read_csv(
    f"ukbiobank/gwas/saige/genetic_similarities/colocalisations/coloc_results_susie.csv")
susie_processed = susie_probable_variant(susie_results=susie_results,
                                         coloc_dict=susie_find_colocs(susie_summary=susie_summaries))

""" Identify genes 'responsible' for colocalisation """


def read_in_gene_set() -> pd.DataFrame:
    """
    Read in Ensembl's hg19 gene set.
    :return:
    """
    ensembl_hg19 = pd.read_csv('https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz', comment='#', sep='\t',
                               header=None, names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                               usecols=['chr', 'feature', 'start', 'end', 'attribute'])[['chr', 'feature', 'start', 'end', 'attribute']]
    ensembl_hg19 = ensembl_hg19[(~ensembl_hg19['chr'].isin(['Y', 'MT'])) & (ensembl_hg19['feature'] == 'gene')].drop(columns='feature').copy()

    # process attribute column to generate multiple columns for each attribute listed
    regex_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    attributes = ensembl_hg19['attribute'].map(
        lambda s: dict(regex_pattern.findall(str(s).strip().strip('"').replace('""','"')))
    ).apply(pd.Series)
    ensembl_hg19 = pd.concat([ensembl_hg19, attributes], axis=1).drop(columns='attribute')

    return ensembl_hg19


def map_snp_to_gene(snp_id: str, gene_set: pd.DataFrame) -> pd.Series:
    """
    Given an SNP ID of form 'chr:pos:ref:alt', identify the Ensembl hg19 gene to which it maps. Mapping means directly
    within the gene or within 2kb of its start or end.
    :param snp_id:
    :param gene_set:
    :return:
    """
    chromo = snp_id.split(':')[0]
    pos = int(snp_id.split(':')[1])

    mapped_genes = gene_set[(gene_set['chr'].astype(str) == chromo) & (gene_set['start'] <= pos + 2000) & (gene_set['end'] >= pos - 2000)].copy()
    mapped_genes_dict = {'ensg_id': '|'.join(mapped_genes.gene_id.to_list()),
                         'gene_symbol': '|'.join(mapped_genes.gene_name.to_list())}

    return pd.Series(mapped_genes_dict)


def map_coloc_results_to_genes(coloc_results_df: pd.DataFrame, gene_set: pd.DataFrame) -> pd.DataFrame:
    """
    Given a set of coloc results which contain SNP IDs in the 'snp' column, map the SNP IDs to Ensembl hg19 genes.
    :param coloc_results_df:
    :param gene_set:
    :return:
    """
    coloc_results_df_working = coloc_results_df.copy()
    coloc_results_df_working[['mapped_ensg_ids', 'mapped_gene_symbols']] = coloc_results_df_working['snp'].apply(
        map_snp_to_gene, gene_set=gene_set)
    coloc_results_df_working['mapped_ensg_ids'] = coloc_results_df_working['mapped_ensg_ids'].replace('', 'Unmapped')
    coloc_results_df_working['mapped_gene_symbols'] = coloc_results_df_working['mapped_gene_symbols'].replace('', 'Unmapped')

    return coloc_results_df_working


ensembl_hg19 = read_in_gene_set()
abf_processed = map_coloc_results_to_genes(coloc_results_df=abf_processed, gene_set=ensembl_hg19)
susie_processed = map_coloc_results_to_genes(coloc_results_df=susie_processed, gene_set=ensembl_hg19)

coloc_results = pd.concat([abf_processed.assign(algo='abf'),
                           susie_processed.assign(algo='susie')])
coloc_results.to_csv(f"~/data/internal/genomics/colocalisations/coloc_abf_susie.csv", index=False)

# then plot a scatterplot-heatmap with ABF top and SuSiE bottom
""" Plot """


def prep_colocs_for_plotting(indiv_coloc_result: pd.DataFrame) -> pd.DataFrame:
    """
    If a disease pair has multiple coloc results, count how many of the 'best'/'strongest' result exist for that pair.
    Also select relevant columns and split 'combo' into separate diseases.
    :param indiv_coloc_result:
    :return:
    """
    indiv_coloc_result_working = indiv_coloc_result[['combo', 'category', 'mapped_gene_symbols']].copy()

    combo_counts = indiv_coloc_result_working.combo.value_counts()
    combos_to_resolve = combo_counts[combo_counts > 1].index.tolist()

    indiv_coloc_result_singlets = indiv_coloc_result_working[
        ~indiv_coloc_result_working['combo'].isin(combos_to_resolve)].copy()
    indiv_coloc_result_singlets['freq'] = 1

    resolved_list = []
    for combo in combos_to_resolve:
        combo_result = indiv_coloc_result_working[indiv_coloc_result_working['combo'] == combo].copy()
        combo_result['category'] = combo_result['category'].astype('category')
        combo_result['category'] = combo_result['category'].cat.set_categories(['strong_coloc-strong_variant',
                                                                                'strong_coloc-moderate_variant',
                                                                                'moderate_coloc-strong_variant',
                                                                                'moderate_coloc-moderate_variant'])
        sorted_categories = combo_result['category'].value_counts().sort_index()
        sorted_categories = sorted_categories[sorted_categories > 0]
        winning_category = sorted_categories.index.tolist()[0]
        joined_mapped_genes = '|'.join(list(set(combo_result[combo_result['category'] == winning_category]['mapped_gene_symbols'])))
        resolved_list.append({'combo': combo,
                              'category': winning_category,
                              'freq': sorted_categories.values[0],
                              'mapped_gene_symbols': joined_mapped_genes})

    processed_data = pd.concat([indiv_coloc_result_singlets, pd.DataFrame(resolved_list)])
    processed_data[['disease1', 'disease2']] = processed_data['combo'].str.split('-', expand=True)

    return processed_data.drop(columns='combo')


disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

plot_path = f"~/ch2_genomics/2.13 colocalisations/plots/colocalisations_abf_susie"

data_to_plot = pd.concat([prep_colocs_for_plotting(abf_processed),
                          prep_colocs_for_plotting(susie_processed).rename(
                              columns={'disease1': 'disease2', 'disease2': 'disease1'})])
data_to_plot = (data_to_plot
                .merge(disease_info[['icd10_three_letter', 'code_chapter']]
                       .rename(columns={'icd10_three_letter': 'disease1'}))
                .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                .merge(disease_info[['icd10_three_letter', 'code_chapter']]
                       .rename(columns={'icd10_three_letter': 'disease2'}))
                .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                          categories=sorted(data_to_plot['disease1'].unique()))
data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                          categories=sorted(data_to_plot['disease2'].unique()))
data_to_plot['category'] = data_to_plot['category'].map({
    'strong_coloc-strong_variant': 'PP.H4.abf > 0.9; \nSNP.PP.H4 > 0.9',
    'strong_coloc-moderate_variant': 'PP.H4.abf > 0.9; \n0.9 >= SNP.PP.H4 > 0.5',
    'moderate_coloc-strong_variant': '0.9 >= PP.H4.abf > 0.5; \nSNP.PP.H4 > 0.9',
    'moderate_coloc-moderate_variant': '0.9 >= PP.H4.abf > 0.5; \n0.9 >= SNP.PP.H4 > 0.5'
})

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

marker_dict = {k: v for k, v in {
        'PP.H4.abf > 0.9; \nSNP.PP.H4 > 0.9': 'o',
        'PP.H4.abf > 0.9; \n0.9 >= SNP.PP.H4 > 0.5': 'v',
        '0.9 >= PP.H4.abf > 0.5; \nSNP.PP.H4 > 0.9': '^',
        '0.9 >= PP.H4.abf > 0.5; \n0.9 >= SNP.PP.H4 > 0.5': 's'
    }.items() if k in data_to_plot.category.unique()}
style_order = [x for x in [
        'PP.H4.abf > 0.9; \nSNP.PP.H4 > 0.9',
        'PP.H4.abf > 0.9; \n0.9 >= SNP.PP.H4 > 0.5',
        '0.9 >= PP.H4.abf > 0.5; \nSNP.PP.H4 > 0.9',
        '0.9 >= PP.H4.abf > 0.5; \n0.9 >= SNP.PP.H4 > 0.5'
    ] if x in data_to_plot.category.unique()]

fig, ax = plt.subplots(figsize=(7, 5))
sns.scatterplot(data=data_to_plot.rename(
    columns={'category': 'Category', 'freq': 'Frequency', 'mapped_gene_symbols': 'Gene'}
),
    x='disease1_as_num', y='disease2_as_num', size='Frequency', style='Category', zorder=3, ax=ax, markers=marker_dict,
    style_order=style_order, sizes=(50, 120), hue='Gene')
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', fontsize=10)
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
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
    if handle._label in ['Gene', 'Frequency', 'Category']:  # section header
        text.set_fontweight('bold')

plt.tight_layout()
plt.savefig(f"{plot_path}.png")
plt.savefig(f"{plot_path}.svg")
plt.close()
