"""
Script to process SAIGE results on both genotyped and imputed data, for either minimal or lifestyles-adjusted results.
"""
from pandas import read_csv, concat, merge, DataFrame
from scipy.stats import chi2

import argparse

parser = argparse.ArgumentParser(description='Reading in which disease to plot.')
parser.add_argument('--adjustment', type=str, required=True, choices=['minimal', 'lifestyles'], help='GWAS adjustment')
args = parser.parse_args()  # parse arguments
adjustment = args.adjustment
print(f"Adjustment: {args.adjustment}", flush=True)

# svat suffix
if adjustment == 'minimal':
    svat_suffix = ''
elif adjustment == 'lifestyles':
    svat_suffix = '_lifestyles'
else:
    raise ValueError('adjustment should be one of minimal or lifestyles.')

# read in disease info
disease_fields_icd10 = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')

# p-val thresholds
pval_threshold = (5 * (10 ** -8)) / 68  # bonferroni corrected
lenient_val_threshold = 5 * (10 ** -8)  # normal threshold

# INFO score thresholds
info_threshold = 0.8
lenient_info_threshold = 0.5

# genotyped data REF/ALT data
geno_ref_alt_info = read_csv('~/data/external/genomics/ukb_geno_ref_alt_info.csv')
# geno_ref_alt_info = read_csv('ukbiobank/gwas/genotype_data/ukb_geno_ref_alt_info.csv')
geno_ref_alt_info.loc[geno_ref_alt_info['chr'] == 23, 'chr'] = 'X'
geno_ref_alt_info['uniqueID_ref_alt'] = geno_ref_alt_info['chr'].astype(str) + '_' + geno_ref_alt_info['pos'].astype(str) + '_' + geno_ref_alt_info['ref'] + '_' + geno_ref_alt_info['alt']
geno_ref_alt_info['ref_alt_sorted'] = geno_ref_alt_info['ref'] + '_' + geno_ref_alt_info['alt']
geno_ref_alt_info['ref_alt_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                       geno_ref_alt_info['ref_alt_sorted'].tolist()]
geno_ref_alt_info['uniqueID_sorted'] = geno_ref_alt_info['chr'].astype(str) + '_' + geno_ref_alt_info['pos'].astype(str) + '_' + geno_ref_alt_info['ref_alt_sorted']

# imputed data REF/ALT data - set up as dictionary. MFI files are ref-first, i.e. A1 = REF, A2 = ALT
imputed_ref_alt = {}
chr_list = list(range(1, 23)) + ['XY']
# chr_list = [21, 'XY']
for i in chr_list:
    mfi_data = read_csv('~/data/external/genomics/ukb_imputed_data/c' + str(i) + '_info_passed.mfi.txt',
                        delim_whitespace=True)
    # mfi_data = read_csv('ukbiobank/gwas/imputed_data/mfi_files/c' + str(i) + '_info_passed.mfi.txt',
    #                     delim_whitespace=True)
    mfi_data['uniqueID_ref_alt'] = str(i) + '_' + mfi_data['pos'].astype(str) + '_' + mfi_data['A1'] + '_' + mfi_data['A2']
    mfi_data['ref_alt_sorted'] = mfi_data['A1'] + '_' + mfi_data['A2']
    mfi_data['ref_alt_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                  mfi_data['ref_alt_sorted'].tolist()]
    mfi_data['uniqueID_sorted'] = str(i) + '_' + mfi_data['pos'].astype(str) + '_' + mfi_data['ref_alt_sorted']
    imputed_ref_alt[str(i)] = mfi_data


""" For each disease, process the per-chromosome results (REF/ALT allele check, p-val check, add data for KG) 
and collate. """
geno_results_list = []
imputed_results_list = []
all_diseases_all_processed_results = []
adjusted_assocs_list = []
for disease in disease_fields_icd10['disease_field'].tolist():
    # disease = 'p130708'
    print('Starting analysis on ' + disease, flush=True)
    icd10_code = disease_fields_icd10[disease_fields_icd10['disease_field'] == disease]['icd10_three_letter'].tolist()[
        0]
    results_dir = '~/ch2_genomics/2.5 full_gwas_runs/saige_sparse_svat' + svat_suffix + '/' + disease + '/'

    disease_all_results = DataFrame()  # collect all unprocessed results for HDL intercept-based correction of chi-squared stats
    """ Process geno results """
    plink_results = read_csv(results_dir + disease + '_pan-ukbb-eur.geno.svat',
                             delim_whitespace=True, dtype={'CHR': str})

    # check whether Allele2 is ALT (0/1/-1 for no/yes/alleles in allele pairs don't match)
    plink_results['alleles_sorted'] = plink_results['Allele1'] + '_' + plink_results['Allele2']
    plink_results['alleles_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                       plink_results['alleles_sorted'].tolist()]
    plink_results['uniqueID_sorted'] = plink_results['CHR'].astype(str) + '_' + plink_results['POS'].astype(
        str) + '_' + plink_results['alleles_sorted']
    plink_results = merge(plink_results, geno_ref_alt_info[['ref', 'alt', 'uniqueID_sorted']],
                          on='uniqueID_sorted', how='inner')
    plink_results['Allele2_is_UKB_alt'] = -1
    plink_results.loc[(plink_results['Allele1'] == plink_results['ref']) &
                      (plink_results['Allele2'] == plink_results['alt']), 'Allele2_is_UKB_alt'] = 1
    plink_results.loc[(plink_results['Allele1'] == plink_results['alt']) &
                      (plink_results['Allele2'] == plink_results['ref']), 'Allele2_is_UKB_alt'] = 0
    plink_results.drop(['alleles_sorted', 'uniqueID_sorted'], axis=1, inplace=True)
    plink_results.rename(columns={'ref': 'UKB_ref', 'alt': 'UKB_alt'}, inplace=True)

    # make PLINK results have the same columns as the imputed results, and add unprocessed results to collation
    plink_results_with_info = plink_results.copy()
    plink_results_with_info['imputationInfo'] = 2.0
    plink_results_with_info['variant_data_type'] = 'genotyped'
    plink_results_with_info.drop('MissingRate', axis=1, inplace=True)
    disease_all_results = concat([disease_all_results, plink_results_with_info])

    # add labels to stratify by p-value and INFO score
    plink_results = plink_results[plink_results['p.value'] < lenient_val_threshold].copy()
    plink_results['assoc_label'] = 'less_assoc'
    plink_results.loc[plink_results['p.value'] < pval_threshold, 'assoc_label'] = 'assoc'

    # add extra required info for KG
    plink_results['method'] = 'SAIGE'
    plink_results['variant_data_type'] = 'genotyped'
    plink_results['source_group'] = 'pan_ukbb_eur_UKB'
    plink_results['source_db'] = 'UKB'
    plink_results['disease'] = icd10_code
    plink_results['assembly'] = 'hg19'
    plink_results['adjustment'] = adjustment

    # filter out MHC co-ords
    plink_results = plink_results[(plink_results['CHR'] != '6') |
                                  (plink_results['POS'] < 28477797) |
                                  (plink_results['POS'] > 33448354)].copy()

    all_chr_imputed_results_list = []
    # chr_list = [21, 'XY']
    for i in chr_list:
        # i = '15'
        print('Starting analysis on chr' + str(i), flush=True)

        """ Process imputed results """
        imputed_results = read_csv(results_dir + disease + '_c' + str(i) + '_pan-ukbb-eur.imputed.svat',
                                   delim_whitespace=True, dtype={'CHR': str})

        # check whether Allele2 is ALT (0/1/-1 for no/yes/alleles in allele pairs don't match)
        imputed_results['alleles_sorted'] = imputed_results['Allele1'] + '_' + imputed_results['Allele2']
        imputed_results['alleles_sorted'] = ['_'.join(sorted(pair.split('_'))) for pair in
                                             imputed_results['alleles_sorted'].tolist()]
        imputed_results['uniqueID_sorted'] = imputed_results['CHR'].astype(str) + '_' + imputed_results[
            'POS'].astype(
            str) + '_' + imputed_results['alleles_sorted']
        imputed_results = merge(imputed_results, imputed_ref_alt[str(i)][['A1', 'A2', 'uniqueID_sorted']],
                                on='uniqueID_sorted', how='inner')
        imputed_results['Allele2_is_UKB_alt'] = -1
        imputed_results.loc[(imputed_results['Allele1'] == imputed_results['A1']) &
                            (imputed_results['Allele2'] == imputed_results['A2']), 'Allele2_is_UKB_alt'] = 1
        imputed_results.loc[(imputed_results['Allele1'] == imputed_results['A2']) &
                            (imputed_results['Allele2'] == imputed_results['A1']), 'Allele2_is_UKB_alt'] = 0
        imputed_results.drop(['alleles_sorted', 'uniqueID_sorted'], axis=1, inplace=True)
        imputed_results.rename(columns={'A1': 'UKB_ref', 'A2': 'UKB_alt'}, inplace=True)

        # add unprocessed results to collation
        imputed_results_with_data_type = imputed_results.copy()
        imputed_results_with_data_type['variant_data_type'] = 'imputed'
        disease_all_results = concat([disease_all_results, imputed_results])

        # filter by MAF - done this way to not assume that Allele2 is the minor allele, i.e. we account for if it's the major allele
        imputed_results = imputed_results[abs(imputed_results['AF_Allele2'] - 0.5) < 0.49].copy()

        # add labels to stratify by p-value and INFO score
        imputed_results = imputed_results[(imputed_results['p.value'] < lenient_val_threshold) &
                                          (imputed_results['imputationInfo'] > lenient_info_threshold)].copy()
        imputed_results['assoc_label'] = 'less_assoc'
        imputed_results.loc[imputed_results['p.value'] < pval_threshold, 'assoc_label'] = 'assoc'
        imputed_results['info_label'] = 'med_info'
        imputed_results.loc[imputed_results['imputationInfo'] > info_threshold, 'info_label'] = 'high_info'

        # add extra required info for KG
        imputed_results['method'] = 'SAIGE'
        imputed_results['variant_data_type'] = 'imputed'
        imputed_results['source_group'] = 'pan_ukbb_eur_UKB'
        imputed_results['source_db'] = 'UKB'
        imputed_results['disease'] = icd10_code
        imputed_results['assembly'] = 'hg19'
        imputed_results['adjustment'] = adjustment

        # filter out MHC
        imputed_results = imputed_results[(imputed_results['CHR'] != '6') |
                                          (imputed_results['POS'] < 28477797) |
                                          (imputed_results['POS'] > 33448354)].copy()

        all_chr_imputed_results_list.append(imputed_results)

    # make per-disease results
    all_chr_imputed_results = concat(all_chr_imputed_results_list)

    plink_results.to_csv(results_dir + disease + '_pan-ukbb-eur_assoc_geno_saige_' + adjustment + '.csv', index=False)
    all_chr_imputed_results.to_csv(results_dir + disease + '_pan-ukbb-eur_assoc_imputed_saige_' + adjustment + '.csv', index=False)

    geno_results_list.append(plink_results)
    imputed_results_list.append(all_chr_imputed_results)

    """ Combine uncorrected geno and imputed processed results into one dataset """
    processed_plink_results_with_info = plink_results.copy()
    processed_plink_results_with_info['imputationInfo'] = 2.0
    processed_plink_results_with_info['info_label'] = 'high_info'
    processed_plink_results_with_info.drop('MissingRate', axis=1, inplace=True)
    # filter dups by UKB allele pair rather than Allele1/2 pair since allele pairs are sometimes duplicated between
    #  genotyped and imputed sets, but with Allele1/2 flipped - want to remove the imputed duplicates
    disease_all_processed_results = concat([processed_plink_results_with_info, all_chr_imputed_results]).drop_duplicates(['CHR', 'POS', 'UKB_ref', 'UKB_alt'])
    all_diseases_all_processed_results.append(disease_all_processed_results)

    """ Remove duplicate alleles from unprocessed results, correct using HDL intercept estimates, and process """
    # remove duplicate alleles
    disease_all_results = disease_all_results.drop_duplicates(['CHR', 'POS', 'UKB_ref', 'UKB_alt']).reset_index(drop=True)

    # convert p-values to chi-squared, correct chi-squareds, convert back to p-values
    chi2_stats = chi2.ppf(1 - disease_all_results['p.value'], df=1)

    # process new adjusted results
    disease_all_results_assoc_info = disease_all_results[(disease_all_results['adjusted_pval'] < pval_threshold) &
                                                         (disease_all_results['imputationInfo'] > info_threshold)].copy()
    disease_all_results_assoc_info['assoc_label'] = 'assoc'
    disease_all_results_assoc_info['info_label'] = 'high_info'
    disease_all_results_assoc_info['method'] = 'SAIGE (HDL intercept-adjusted)'
    disease_all_results_assoc_info['variant_data_type'] = 'imputed'
    disease_all_results_assoc_info.loc[disease_all_results_assoc_info['imputationInfo'] == 2, 'variant_data_type'] = 'genotyped'
    disease_all_results_assoc_info['source_group'] = 'pan_ukbb_eur_UKB'
    disease_all_results_assoc_info['source_db'] = 'UKB'
    disease_all_results_assoc_info['disease'] = icd10_code
    disease_all_results_assoc_info['assembly'] = 'hg19'
    disease_all_results_assoc_info['adjustment'] = adjustment
    disease_all_results_assoc_info = disease_all_results_assoc_info[(disease_all_results_assoc_info['CHR'] != '6') |
                                                                    (disease_all_results_assoc_info['POS'] < 28477797) |
                                                                    (disease_all_results_assoc_info['POS'] > 33448354)].copy()
    adjusted_assocs_list.append(disease_all_results_assoc_info)

# collate into pan-disease pan-chr results
geno_results = concat(geno_results_list)
imputed_results = concat(imputed_results_list)

geno_results.to_csv('~/ch2_genomics/2.5 full_gwas_runs/pan-ukbb-eur_assoc_geno_saige_' + adjustment + '.csv', index=False)
imputed_results.to_csv('~/ch2_genomics/2.5 full_gwas_runs/pan-ukbb-eur_assoc_imputed_saige_' + adjustment + '.csv', index=False)

# export the combined geno-imputed results
concat(all_diseases_all_processed_results).to_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_' + adjustment + '.csv', index=False)
