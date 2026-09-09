"""
Script to process REGENIE results on both genotyped and imputed data.
"""
from pandas import read_csv, concat, merge
from math import log10

# read in disease info
disease_fields = read_csv('~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_disease_fields.txt')
disease_fields_icd10 = read_csv('~/data/internal/phenotype_coding/disease_fields_icd10_info.csv')

# p-val thresholds
pval_threshold = 5 * (10 ** -8)
lenient_pval_threshold = 5 * (10 ** -7)

# -log10 p-val thresholds
log_pval_threshold = log10(pval_threshold) * -1
log_lenient_pval_threshold = log10(lenient_pval_threshold) * -1

# MAF threshold
maf_threshold = 0.01

# genotyped data REF/ALT data
geno_ref_alt_info = read_csv('~/data/external/genomics/ukb_geno_ref_alt_info.csv')
geno_ref_alt_info_for_allele_check = geno_ref_alt_info.assign(
    pair=geno_ref_alt_info.apply(lambda x: f"{x['ref']}/{x['alt']}", axis=1))
geno_ref_alt_info_for_allele_check['pair'] = ['/'.join(sorted(pair.split('/'))) for pair in
                                              geno_ref_alt_info_for_allele_check['pair'].tolist()]

# imputed data REF/ALT data - set up as dictionary
imputed_ref_alt = {}
chr_list = list(range(1, 23)) + ['XY']
# chr_list = [19]
for i in chr_list:
    mfi_data = read_csv('~/data/external/genomics/ukb_imputed_data/c' + str(i) + '_info_passed.mfi.txt',
                        low_memory=False,
                        delim_whitespace=True)
    mfi_data_for_allele_check = mfi_data.assign(
        pair=mfi_data.apply(lambda x: f"{x['A1']}/{x['A2']}", axis=1))
    mfi_data_for_allele_check['pair'] = ['/'.join(sorted(pair.split('/'))) for pair in
                                         mfi_data_for_allele_check['pair'].tolist()]
    imputed_ref_alt[str(i)] = mfi_data_for_allele_check

""" For each disease, process the genotyped and imputed per-chromosome results (MAF filter, REF/ALT allele check, 
p-val check, add data for KG) and collate. """
geno_results_list = []
imputed_results_list = []
for disease in disease_fields['disease_field'].tolist():
    # disease = 'p131036'
    icd10_code = disease_fields_icd10[disease_fields_icd10['disease_field'] == disease]['icd10_three_letter'].tolist()[
        0]
    results_dir = 'ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step2/' + disease + '/'

    """ Process geno results """
    geno_results = read_csv(results_dir + disease + '_pan-ukbb-eur.geno.step2_' + disease + '.regenie', delim_whitespace=True)

    # filter by MAF > 0.01
    geno_results_maf = geno_results[(geno_results['A1FREQ'] > maf_threshold) &
                                    (geno_results['A1FREQ'] < 1 - maf_threshold)].copy()

    # add labels to stratify by p-value
    geno_results_maf_assoc = geno_results_maf[geno_results_maf['LOG10P'] > log_pval_threshold].copy()
    geno_results_maf_assoc['assoc_label'] = 'assoc'
    geno_results_maf_less_assoc = geno_results_maf[(geno_results_maf['LOG10P'] <= log_pval_threshold) &
                                                   (geno_results_maf['LOG10P'] > log_lenient_pval_threshold)].copy()
    geno_results_maf_less_assoc['assoc_label'] = 'less_assoc'
    geno_results_maf_assoc_label = concat([geno_results_maf_assoc, geno_results_maf_less_assoc])

    # add extra required info for KG
    geno_results_maf_assoc_label['method'] = 'REGENIE'
    geno_results_maf_assoc_label['variant_data_type'] = 'genotyped'
    geno_results_maf_assoc_label['source_group'] = 'pan_ukbb_eur_UKB'
    geno_results_maf_assoc_label['source_db'] = 'UKB'
    geno_results_maf_assoc_label['disease'] = icd10_code
    geno_results_maf_assoc_label['assembly'] = 'hg19'

    # filter out MHC
    geno_results_maf_assoc_label_no_mhc = geno_results_maf_assoc_label[(geno_results_maf_assoc_label['CHROM'] != 6) |
                                                                       (geno_results_maf_assoc_label['GENPOS'] < 28477797) |
                                                                       (geno_results_maf_assoc_label['GENPOS'] > 33448354)]
    # REF/ALT allele check
    geno_results_for_allele_check = geno_results_maf_assoc_label_no_mhc.assign(
        gwas_pair=geno_results_maf_assoc_label_no_mhc.apply(lambda x: f"{x['ALLELE0']}/{x['ALLELE1']}", axis=1))
    geno_results_for_allele_check['gwas_pair'] = ['/'.join(sorted(pair.split('/'))) for pair in
                                                 geno_results_for_allele_check['gwas_pair'].tolist()]
    geno_results_ref_alt = merge(geno_results_for_allele_check, geno_ref_alt_info_for_allele_check,
                                 left_on=['ID', 'gwas_pair'],
                                 right_on=['id', 'pair'])
    geno_results_allele_checked = geno_results_maf_assoc_label_no_mhc.copy()
    geno_results_allele_checked['ALLELE0'] = geno_results_ref_alt['ref'].tolist()
    geno_results_allele_checked['ALLELE1'] = geno_results_ref_alt['alt'].tolist()
    geno_results_allele_checked.rename(columns={'ALLELE0': 'REF',
                                                'ALLELE1': 'ALT'},
                                       inplace=True)
    geno_results_allele_checked.drop(['A1FREQ'], axis=1, inplace=True)

    """ Process imputed results """
    all_chr_imputed_results_list = []
    for i in chr_list:
        # i = '19'

        imputed_raw_results = read_csv(results_dir + disease + '_c' + str(i) + '_pan-ukbb-eur.imputed.step2_' + disease + '.regenie',
                                       low_memory=False,
                                       delim_whitespace=True)

        imputed_maf = imputed_raw_results[(imputed_raw_results['A1FREQ'] > maf_threshold) &
                                          (imputed_raw_results['A1FREQ'] < 1 - maf_threshold)].copy()

        # add labels to stratify by p-value and INFO score
        imputed_maf_assoc = imputed_maf[imputed_maf['LOG10P'] > log_pval_threshold].copy()
        imputed_maf_assoc['assoc_label'] = 'assoc'
        imputed_maf_less_assoc = imputed_maf[(imputed_maf['LOG10P'] > log_lenient_pval_threshold) &
                                             (imputed_maf['LOG10P'] <= log_pval_threshold)].copy()
        imputed_maf_less_assoc['assoc_label'] = 'less_assoc'
        imputed_maf_assoc_label = concat([imputed_maf_assoc, imputed_maf_less_assoc])

        imputed_maf_assoc_label_med_info = imputed_maf_assoc_label[(imputed_maf_assoc_label['INFO'] > 0.5) &
                                                                   (imputed_maf_assoc_label['INFO'] <= 0.8)].copy()
        imputed_maf_assoc_label_med_info['info_label'] = 'med_info'
        imputed_maf_assoc_label_high_info = imputed_maf_assoc_label[imputed_maf_assoc_label['INFO'] > 0.8].copy()
        imputed_maf_assoc_label_high_info['info_label'] = 'high_info'
        imputed_maf_assoc_info_labels = concat([imputed_maf_assoc_label_high_info, imputed_maf_assoc_label_med_info])

        # add extra required info for KG
        imputed_maf_assoc_info_labels['method'] = 'REGENIE'
        imputed_maf_assoc_info_labels['variant_data_type'] = 'imputed'
        imputed_maf_assoc_info_labels['source_group'] = 'pan_ukbb_eur_UKB'
        imputed_maf_assoc_info_labels['source_db'] = 'UKB'
        imputed_maf_assoc_info_labels['disease'] = icd10_code
        imputed_maf_assoc_info_labels['assembly'] = 'hg19'

        # filter out MHC
        imputed_maf_assoc_info_labels_no_mhc = imputed_maf_assoc_info_labels[(imputed_maf_assoc_info_labels['CHROM'] != 6) |
                                                                             (imputed_maf_assoc_info_labels['GENPOS'] < 28477797) |
                                                                             (imputed_maf_assoc_info_labels['GENPOS'] > 33448354)].copy()

        # REF/ALT allele check
        imputed_results_for_allele_check = imputed_maf_assoc_info_labels_no_mhc.assign(
            pair=imputed_maf_assoc_info_labels_no_mhc.apply(lambda x: f"{x['ALLELE0']}/{x['ALLELE1']}", axis=1))
        imputed_results_for_allele_check['pair'] = ['/'.join(sorted(pair.split('/'))) for pair in
                                                    imputed_results_for_allele_check['pair'].tolist()]
        imputed_for_ref_alt = merge(imputed_results_for_allele_check, imputed_ref_alt[str(i)],
                                    left_on=['ID', 'pair'],
                                    right_on=['snpID', 'pair'])
        imputed_results_allele_checked = imputed_maf_assoc_info_labels_no_mhc.copy()
        imputed_results_allele_checked['ALLELE0'] = imputed_for_ref_alt['A1'].tolist()
        imputed_results_allele_checked['ALLELE1'] = imputed_for_ref_alt['A2'].tolist()
        imputed_results_allele_checked.rename(columns={'ALLELE0': 'REF',
                                                       'ALLELE1': 'ALT'}, inplace=True)
        imputed_results_allele_checked.drop(['A1FREQ'], axis=1, inplace=True)
        all_chr_imputed_results_list.append(imputed_results_allele_checked.sort_values(by=['CHROM', 'GENPOS']))

    # make per-disease results
    all_chr_geno_results = geno_results_allele_checked.copy()
    all_chr_imputed_results = concat(all_chr_imputed_results_list)

    all_chr_geno_results.to_csv(results_dir + disease + '_pan-ukbb-eur_assoc_geno_regenie.csv', index=False)
    all_chr_imputed_results.to_csv(results_dir + disease + '_pan-ukbb-eur_assoc_imputed_regenie.csv', index=False)

    geno_results_list.append(all_chr_geno_results)
    imputed_results_list.append(all_chr_imputed_results)

# collate into pan-disease pan-chr results
geno_results = concat(geno_results_list)
imputed_results = concat(imputed_results_list)

geno_results.to_csv('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/pan-ukbb-eur_assoc_geno_regenie.csv', index=False)
imputed_results.to_csv('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/pan-ukbb-eur_assoc_imputed_regenie.csv', index=False)
