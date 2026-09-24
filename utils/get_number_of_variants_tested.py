""" Get total number of variants tested in the GWAS, per-disease and per-disease pair. """
import pandas as pd


def remove_variants_from_genotyped_variants(disease_field: str, genotyped_variants: pd.DataFrame) -> pd.DataFrame:
    """
    Some variants that should be removed from the genotyped set are disease-specific, and are removed here.
    :param disease_field:
    :param genotyped_variants:
    :return:
    """
    variants_to_remove = pd.read_csv(f"~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/variants_to_exclude/{disease_field}_autosomes-cX-vars_to_exclude.tsv", sep='\t', header=None)
    genotyped_variants_post_removal = genotyped_variants[
        (~genotyped_variants['id'].isin(variants_to_remove[0].to_list())) &
        (~genotyped_variants['chr'].isin(['24', '25', '26']))]

    return genotyped_variants_post_removal


def read_in_and_process_information(logger):
    """
    Read in data.
    :return: Dict of (1) disease info
    (2) dictionary of pd.Series of unique genotyped variant IDs per disease
    (3) pd.Series of unique imputed variant IDs.
    """
    """ Read in disease info """
    disease_info = pd.read_csv(
        '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    field_to_code_mapper = dict(zip(disease_info['disease_field'], disease_info['icd10_three_letter']))

    """ Read in genotyped variants and get tested set for each disease """
    genotyped_info = pd.read_csv('~/data/external/genomics/ukb_geno_ref_alt_info.csv', dtype={'chr': str})
    genotyped_info = genotyped_info[~genotyped_info['chr'].isin(['25', '26'])].copy()
    genotyped_info.loc[genotyped_info['chr'] == '23', 'chr'] = 'X'
    genotyped_info.loc[genotyped_info['chr'] == '24', 'chr'] = 'XY'  # X chromosome PAR
    genotyped_info['alleles_sorted'] = genotyped_info['ref'] + '_' + genotyped_info['alt']
    genotyped_info['alleles_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                        genotyped_info['alleles_sorted'].tolist()]
    genotyped_info['uniqueID_sorted'] = genotyped_info['chr'].astype(str) + '_' + genotyped_info['pos'].astype(
        str) + '_' + genotyped_info['alleles_sorted']

    genotyped_dict = {}
    for disease in disease_info.disease_field:
        disease_genotyped_info = remove_variants_from_genotyped_variants(disease_field=disease,
                                                                         genotyped_variants=genotyped_info)
        genotyped_dict[field_to_code_mapper[disease]] = disease_genotyped_info['uniqueID_sorted']
        logger.info(f"Read in genotyped variants for {disease}: {len(disease_genotyped_info['uniqueID_sorted'])} variants.")

    """ For every chromosome, read in the imputed data, isolate the relevant info """

    n_vars_total = 0
    imputed_unique_ids = pd.Series()
    chromo_list = sorted(list(genotyped_info.chr.unique()))
    logger.info(f"Reading in imputed variants for {', '.join(['c' + x for x in chromo_list])}")
    for chromo in chromo_list:
        chromo_imputed_info = pd.read_csv(f"~/data/external/genomics/ukb_imputed_data/mfi_files/c{chromo}_info_passed.mfi.txt",
                                          sep='\s+', usecols=['MAF', 'INFO', 'A1', 'A2', 'pos'])
        # filter on INFO > 0.8 and MAF > 0.01
        chromo_imputed_info = chromo_imputed_info[
            (abs(chromo_imputed_info['MAF'] - 0.5) < 0.49) & (chromo_imputed_info['INFO'] > 0.8)]
        chromo_imputed_info['ref_alt_sorted'] = chromo_imputed_info['A1'] + '_' + chromo_imputed_info['A2']
        chromo_imputed_info['ref_alt_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                                 chromo_imputed_info['ref_alt_sorted'].tolist()]
        chromo_imputed_info['uniqueID_sorted'] = str(chromo) + '_' + chromo_imputed_info['pos'].astype(str) + '_' + \
                                                 chromo_imputed_info['ref_alt_sorted']
        imputed_unique_ids = pd.concat([imputed_unique_ids, chromo_imputed_info['uniqueID_sorted']])
        logger.info(f"Read in imputed variants for c{chromo}")

    return {'disease_info': disease_info, 'genotype_vars_by_disease': genotyped_dict, 'imputed_vars': imputed_unique_ids}


def get_vars_tested_per_disease(disease: str, data_dict) -> pd.Series:
    """
    Get the set of variants tested for a given disease.
    :param disease:
    :param data_dict:
    :return:
    """

    genotyped_vars = data_dict['genotype_vars_by_disease'][disease]
    return pd.Series(pd.concat([genotyped_vars, data_dict['imputed_vars']]).unique())


def get_vars_tested_per_disease_pair(disease1: str, disease2: str, data_dict) -> pd.Series:
    """
    Get the Union set of variants tested for a given disease pair.
    :param disease1:
    :param disease2:
    :param data_dict:
    :return:
    """

    # disease1
    disease1_vars = get_vars_tested_per_disease(disease=disease1, data_dict=data_dict)
    disease2_vars = get_vars_tested_per_disease(disease=disease2, data_dict=data_dict)
    disease1_disease2_vars = pd.Series(pd.concat([disease1_vars, disease2_vars]).unique())

    return disease1_disease2_vars
