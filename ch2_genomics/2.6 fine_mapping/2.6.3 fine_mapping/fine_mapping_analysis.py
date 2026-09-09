"""
Map fine mapped variants to Ensembl hg19 genes and analyse with that info.
"""
from typing import Dict

import collections
import os
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('TkAgg')

adjustment = 'minimal'


def read_in_gene_set() -> pd.DataFrame:
    """
    Read in Ensembl's hg19 gene set.
    :return:
    """
    ensembl_hg19 = \
        pd.read_csv('https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz',
                    comment='#', sep='\t',
                    header=None,
                    names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                    usecols=['chr', 'feature', 'start', 'end', 'attribute'])[
            ['chr', 'feature', 'start', 'end', 'attribute']]
    ensembl_hg19 = ensembl_hg19[(~ensembl_hg19['chr'].isin(['Y', 'MT'])) & (ensembl_hg19['feature'] == 'gene')].drop(
        columns='feature').copy()

    # process attribute column to generate multiple columns for each attribute listed
    regex_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    attributes = ensembl_hg19['attribute'].map(
        lambda s: dict(regex_pattern.findall(str(s).strip().strip('"').replace('""', '"')))
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

    mapped_genes = gene_set[(gene_set['chr'].astype(str) == chromo) & (gene_set['start'] <= pos + 2000) & (
            gene_set['end'] >= pos - 2000)].copy()
    mapped_genes_dict = {'ensg_id': '|'.join(mapped_genes.gene_id.to_list()),
                         'gene_symbol': '|'.join(mapped_genes.gene_name.to_list())}

    return pd.Series(mapped_genes_dict)


def map_mapping_results_to_genes(mapping_results_df: pd.DataFrame, gene_set: pd.DataFrame) -> pd.DataFrame:
    """
    Given a set of fine mapping results which contain SNP IDs in the 'snp' column, map the SNP IDs to Ensembl hg19 genes.
    :param mapping_results_df:
    :param gene_set:
    :return:
    """
    mapping_results_df_working = mapping_results_df.copy()
    mapping_results_df_working[['mapped_ensg_ids', 'mapped_gene_symbols']] = mapping_results_df_working[
        'variant'].apply(
        map_snp_to_gene, gene_set=gene_set)
    mapping_results_df_working['mapped_ensg_ids'] = mapping_results_df_working['mapped_ensg_ids'].replace('',
                                                                                                          'Unmapped')
    mapping_results_df_working['mapped_gene_symbols'] = mapping_results_df_working['mapped_gene_symbols'].replace('',
                                                                                                                  'Unmapped')

    return mapping_results_df_working


def merge_overlapping_ld_regions(disease_algo_df: pd.DataFrame) -> pd.DataFrame:
    # disease_algo_df = pd.DataFrame for each disease-algo combination, e.g. E11-SuSiE, E11-ABF, etc.
    # aim: map chr:ld_region_start:ld_region_end to a given block

    disease_algo_df = disease_algo_df.sort_values(by=['chr', 'ld_region_start', 'ld_region_end']).reset_index(
        drop=True).copy()

    ld_region_start_store = []
    ld_region_end_store = []
    merged_ld_blocks = list()
    for i in range(len(disease_algo_df)):
        current_chr = disease_algo_df.at[i, 'chr']
        current_start = disease_algo_df.at[i, 'ld_region_start']
        current_end = disease_algo_df.at[i, 'ld_region_end']

        ld_region_start_store.append(current_start)
        ld_region_end_store.append(current_end)

        if i + 1 == len(disease_algo_df):
            # generate the merged ld block and store
            merged_ld_block = '_'.join(
                [str(current_chr), str(min(ld_region_start_store)), str(max(ld_region_end_store))])
            merged_ld_blocks += [merged_ld_block] * len(ld_region_start_store)
            # merged_ld_blocks[min(idx_store):] = merged_ld_block

            ld_region_start_store = []
            ld_region_end_store = []

            continue

        next_chr = disease_algo_df.at[i + 1, 'chr']
        next_start = disease_algo_df.at[i + 1, 'ld_region_start']

        # if there will be a new chromosome in the next step, or there is sufficient gap between regions -> close counting
        if (current_chr != next_chr) or (current_end + 1000000 < next_start):
            # generate the merged ld block and store
            merged_ld_block = '_'.join(
                [str(current_chr), str(min(ld_region_start_store)), str(max(ld_region_end_store))])
            merged_ld_blocks += [merged_ld_block] * len(ld_region_start_store)

            ld_region_start_store = []
            ld_region_end_store = []

    disease_algo_df['merged_ld_region'] = merged_ld_blocks

    return disease_algo_df


def new_credible_set_for_merged_block(merged_block_df: pd.DataFrame) -> pd.DataFrame:
    # for an algo-disease-merged block combination, create a new 'credible set' by reordering variants by pip and dropping duplicates

    merged_block_df2 = merged_block_df.drop(columns='disease_credible_set').sort_values(by='pip',
                                                                                        ascending=False).drop_duplicates(
        subset='variant').reset_index(drop=True).copy()
    merged_block_df2['disease_credible_set'] = [x + 1 for x in list(merged_block_df2.index)]
    merged_block_df2 = merged_block_df2[list(merged_block_df.columns)]

    return merged_block_df2


ensembl_hg19 = read_in_gene_set()
mapped_vars = pd.read_csv(f"~/data/internal/genomics/fine_mapping/{adjustment}/peak_pip0.5_fine_mapped_variants.csv")
mapped_vars = map_mapping_results_to_genes(mapping_results_df=mapped_vars, gene_set=ensembl_hg19)

# to collapse same blocks/genes/variants into a single entry
# for each algo-disease combination, ensure one set of non-overlapping LD regions
mapped_vars_with_merged_regions = pd.concat(
    [merge_overlapping_ld_regions(df) for _, df in mapped_vars.groupby(['disease', 'algo'])])

# for each algo-disease-merged block combination, create new credible set
mapped_vars_with_merged_regions_recredibled = pd.concat([new_credible_set_for_merged_block(df) for _, df in
                                                         mapped_vars_with_merged_regions.groupby(
                                                             ['disease', 'algo', 'merged_ld_region'])])
mapped_vars_with_merged_regions_recredibled.to_csv(
    f"~/data/internal/genomics/fine_mapping/{adjustment}/peak_pip0.5_fine_mapped_variants_blocksmerged_recredibled.csv",
    index=False)

""" High-level statistics """

processed_fine_mapping = pd.read_csv(
    f"~/data/internal/genomics/fine_mapping/{adjustment}/peak_pip0.5_fine_mapped_variants_blocksmerged_recredibled.csv")
processed_fine_mapping['mapped_gene_symbols'] = processed_fine_mapping['mapped_gene_symbols'].str.split('|')
processed_fine_mapping = processed_fine_mapping.explode('mapped_gene_symbols')

susie = processed_fine_mapping[processed_fine_mapping['algo'] == 'SuSiE'].copy()
abf = processed_fine_mapping[processed_fine_mapping['algo'] == 'ABF'].copy()

# get number of diseases
susie[susie.pip > 0.5].disease.nunique()
susie[susie.pip > 0.9].disease.nunique()

# get number of genes and gene frequency
x = (susie[susie.pip > 0.5]
     .assign(mapped_gene_symbols=susie[susie.pip > 0.5]
             .mapped_gene_symbols.str.split('|'))
     .explode('mapped_gene_symbols')
     .drop_duplicates(['disease', 'mapped_gene_symbols']))
y = x.mapped_gene_symbols.value_counts()
yy = x.disease.value_counts()

# number of fine-mapped regions per disease (high confidence fine mapping)
z = pd.DataFrame(susie[susie.pip > 0.9].value_counts(
    ['disease', 'chr', 'merged_ld_region', 'mapped_gene_symbols'])).reset_index().value_counts('disease')

# instances where there are multiple variants fine mapped to a single region
a = susie.value_counts(['disease', 'merged_ld_region'])
