""" Create tables for which genes and variants are significantly shared between ARDs. """
import pandas as pd
import sys
import itertools

from pyprind import ProgBar


def get_entity_overlap(
        data_df: pd.DataFrame, disease_col: str, entity_col: str
) -> pd.DataFrame:
    """
    For each disease pair in disease_pairs, identify the set of entities shared between the two.
    :param data_df:
    :param disease_col:
    :param entity_col:
    :return:
    """
    # shared_entity_dict = collections.defaultdict(list)
    disease_pairs = list(itertools.combinations(sorted(list(data_df[disease_col].unique())), 2))
    shared_entity_dicts = []
    bar = ProgBar(len(disease_pairs), stream=sys.stdout, title='Identifying shared entities')
    for disease_pair in disease_pairs:
        disease1 = disease_pair[0]
        disease2 = disease_pair[1]

        shared_entities = list(set(data_df[data_df[disease_col] == disease1][entity_col].unique()) &
                               set(data_df[data_df[disease_col] == disease2][entity_col].unique()))
        if len(shared_entities) > 0:
            shared_entity_dicts.append({'disease1': disease1, 'disease2': disease2,
                                        'shared_entities': ', '.join(shared_entities)})
        # shared_entity_dict[disease_pair] = shared_entities
        bar.update()

    return pd.DataFrame(shared_entity_dicts)


""" Aux """
disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(str) + ')'

""" Variants """
variants = pd.read_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_lifestyles.csv',
                       dtype={'CHR': str, 'POS': int})
variants = variants[(variants['assoc_label'] == 'assoc') & (variants['info_label'] == 'high_info')]
variants['alleles_sorted'] = variants['UKB_ref'] + '_' + variants['UKB_alt']
variants['alleles_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                              variants['alleles_sorted'].tolist()]
variants['uniqueID_sorted'] = variants['CHR'].astype(str) + '_' + variants['POS'].astype(
    str) + '_' + variants['alleles_sorted']
variants.loc[variants['Allele2_is_UKB_alt'] == 0, 'BETA'] = -1 * variants.loc[
    variants['Allele2_is_UKB_alt'] == 0, 'BETA']
variants = (variants[['disease', 'uniqueID_sorted', 'BETA']]
            .drop_duplicates(subset=['disease', 'uniqueID_sorted'])
            .reset_index(drop=True))
variants_overlap = get_entity_overlap(data_df=variants, disease_col='disease', entity_col='uniqueID_sorted')
# read in overlap significance results
variant_fisher_sig = pd.read_csv('~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_cossims.csv')
variants_overlap_sig = variants_overlap.merge(variant_fisher_sig[variant_fisher_sig['fdr_sig'] == 1][['disease1', 'disease2']])
variants_overlap_sig[['disease1', 'disease2', 'shared_entities']].to_csv('~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_table.csv', index=False)

""" Genes """
all_genes = pd.read_csv('~/data/external/genomics/magma/NCBI37.3/NCBI37.3.gene.loc',
                        sep='\s+', header=None)
genes_assoc = pd.read_csv('~/data/internal/genomics/magma/lifestyles/gene_level_results_assoc.csv')
genes_assoc = genes_assoc[genes_assoc['bonf_sig'] == 1][['GENE', 'ZSTAT', 'disease']]
genes_assoc = genes_assoc.merge(all_genes[[0, 5]].rename(columns={0: 'GENE', 5: 'gene_name'})).drop(columns='GENE').rename(columns={'gene_name': 'GENE'})
genes_overlap = get_entity_overlap(data_df=genes_assoc, disease_col='disease', entity_col='GENE')
# read in overlap significance results
gene_fisher_sig = pd.read_csv('~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher.csv')
genes_overlap_sig = genes_overlap.merge(gene_fisher_sig[gene_fisher_sig['fdr_sig'] == 1][['disease1', 'disease2']])
genes_overlap_sig = (genes_overlap_sig
                .merge(disease_info[['icd10_three_letter', 'disease_field']]
                       .rename(columns={'disease_field': 'disease1'}))
                .drop(columns='disease1').rename(columns={'icd10_three_letter': 'disease1'})
                .merge(disease_info[['icd10_three_letter', 'disease_field']]
                       .rename(columns={'disease_field': 'disease2'}))
                .drop(columns='disease2').rename(columns={'icd10_three_letter': 'disease2'}))
genes_overlap_sig[['disease1', 'disease2', 'shared_entities']].to_csv('~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher_table.csv', index=False)
