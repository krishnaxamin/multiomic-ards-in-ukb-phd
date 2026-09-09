"""
Creates tsvs listing, for each disease and cohort,:
- autosomal variants to exclude (c1-22, cXY)
- cX variants to exclude
- a combined list (autosomal + cX) of SNPs to exclude.
Uses as input stats calculated using PLINK:
- allele frequency (--freq -> .frq)
- variant missingness (--missing -> .lmiss)
- sample missingness (--missing -> .imiss).
(Missingness is calculated from N_MISS and N_GENO in the .lmiss and .imiss files. This is because PLINK calculates
missingness rate to 4 sig figs, so in some instances a variant that is actually above the missingness threshold is
counted as not so because its missingness rate is rounded down to 0.1. Therefore, the missingness rate is calculated
from the raw N_MISS and N_GENO counts - this rate is unaffected by rounding. Reasoning given below at first instance).
Filter thresholds taken from:
- Donertas et al, 2021 (general MAF and HWE)
- Marees et al, 2018 (general lmiss and imiss [also HWE on controls])
- Konig et al, 2014 (application of thresholds for cX variants).

A different HWE threshold is no longer used for cX, since HWE stats are now calculated on cX using PLINK2 and using both
male and female control populations, in accordance with PLINK2's usage.
See this post https://groups.google.com/g/plink2-users/c/vA61ymOtZ5o/m/pBCKizX3AQAJ and 
this paper https://pubmed.ncbi.nlm.nih.gov/27071844/.

From 2024-02-10, calculations are performed sequentially after each filtering step, following the order in Marees et al,
 2018 (Genomic QC 2).
"""
import sys

from pandas import read_csv, concat, DataFrame
from pyprind import ProgBar

# filter thresholds
maf_threshold = 0.01  # filter out maf < 0.01
vmiss_threshold = 0.02  # filter out lmiss > 0.02
smiss_threshold = 0.02  # filter out smiss > 0.02
hwe_autosomes_threshold = 10 ** -6  # filter out hwe < 1e-6
# hwe_cx_threshold = 10 ** -4  # filter out hwe < 1e-4

""" Pipeline after 2024-02-10, after which calculations are performed after each filtering step (Genomic QC 2) """

# 1. Calculate variant missingness

# 2. Remove variants based on missingness
autosomes_vmiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step1_variant_missingness/autosomes_panukbb_eur.vmiss', delim_whitespace=True)
autosomes_vmiss['F_MISS_CALCED'] = (autosomes_vmiss['MISSING_CT'] / autosomes_vmiss['OBS_CT'])
autosomes_vmiss_fail = autosomes_vmiss[autosomes_vmiss['F_MISS_CALCED'] > vmiss_threshold]

cx_males_vmiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step1_variant_missingness/cX_panukbb_eur_male.vmiss', delim_whitespace=True)
cx_males_vmiss['F_MISS_CALCED'] = (cx_males_vmiss['MISSING_CT'] / cx_males_vmiss['OBS_CT'])
cx_females_vmiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step1_variant_missingness/cX_panukbb_eur_female.vmiss', delim_whitespace=True)
cx_females_vmiss['F_MISS_CALCED'] = (cx_females_vmiss['MISSING_CT'] / cx_females_vmiss['OBS_CT'])
cx_mf_diff_vmiss = DataFrame({'ID': cx_females_vmiss['ID'],
                             'male_female_diff': abs(cx_females_vmiss['F_MISS_CALCED'] - cx_males_vmiss['F_MISS_CALCED'])})
cx_vmiss_fail = concat([cx_females_vmiss[cx_females_vmiss['F_MISS_CALCED'] > vmiss_threshold],
                        cx_males_vmiss[cx_males_vmiss['F_MISS_CALCED'] > vmiss_threshold],
                        cx_mf_diff_vmiss[cx_mf_diff_vmiss['male_female_diff'] > vmiss_threshold]])

# export high missingness variants
high_missingness_variants = concat([autosomes_vmiss_fail[['ID']], cx_vmiss_fail[['ID']]]).drop_duplicates(ignore_index=True)
high_missingness_variants.to_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step2_high_missingness_variants_panukbb_eur.tsv', index=False, header=False, sep='\t')

# 3. Calculate sample missingness

# 4. Remove samples based on missingness
autosomes_smiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step3_sample_missingness/autosomes_panukbb_eur.smiss', delim_whitespace=True)
autosomes_smiss['F_MISS_CALCED'] = (autosomes_smiss['MISSING_CT'] / autosomes_smiss['OBS_CT'])
autosomes_smiss_fail = autosomes_smiss[autosomes_smiss['F_MISS_CALCED'] > smiss_threshold]

cx_males_smiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step3_sample_missingness/cX_panukbb_eur_male.smiss', delim_whitespace=True)
cx_males_smiss['F_MISS_CALCED'] = (cx_males_smiss['MISSING_CT'] / cx_males_smiss['OBS_CT'])
cx_females_smiss = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step3_sample_missingness/cX_panukbb_eur_female.smiss', delim_whitespace=True)
cx_females_smiss['F_MISS_CALCED'] = (cx_females_smiss['MISSING_CT'] / cx_females_smiss['OBS_CT'])
cx_smiss_fail = concat([cx_males_smiss[cx_males_smiss['F_MISS_CALCED'] > smiss_threshold],
                        cx_females_smiss[cx_females_smiss['F_MISS_CALCED'] > smiss_threshold]])

high_missingness_samples = concat([autosomes_smiss_fail[['#FID', 'IID']],
                                   cx_smiss_fail[['#FID', 'IID']]]).drop_duplicates(ignore_index=True)
high_missingness_samples['pheno'] = 1
high_missingness_samples.to_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv', index=False, header=False, sep='\t')

""" With the QC-ed sample set, go establish what the ARDs are, and then come back to this point"""

# 5. Calculate MAF

# 6. Remove variants based on MAF
autosomes_freq = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step5_maf/autosomes_panukbb_eur.acount', delim_whitespace=True)
autosomes_freq['ALT_FREQ'] = (autosomes_freq['ALT_CTS'] / autosomes_freq['OBS_CT'])
autosomes_freq_fail = autosomes_freq[abs(autosomes_freq['ALT_FREQ'] - 0.5) > 0.5 - maf_threshold]

cx_males_freq = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step5_maf/cX_panukbb_eur_male.acount', delim_whitespace=True)
cx_males_freq['ALT_FREQ'] = (cx_males_freq['ALT_CTS'] / cx_males_freq['OBS_CT'])
cx_females_freq = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step5_maf/cX_panukbb_eur_female.acount', delim_whitespace=True)
cx_females_freq['ALT_FREQ'] = (cx_females_freq['ALT_CTS'] / cx_females_freq['OBS_CT'])
cx_freq_fail = concat([cx_males_freq[abs(cx_males_freq['ALT_FREQ'] - 0.5) > 0.5 - maf_threshold],
                       cx_females_freq[abs(cx_females_freq['ALT_FREQ'] - 0.5) > 0.5 - maf_threshold]]).drop_duplicates(ignore_index=True)

low_maf_variants = concat([autosomes_freq_fail[['ID']], cx_freq_fail[['ID']]]).drop_duplicates(ignore_index=True)
low_maf_variants.to_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step6_low_maf_variants_panukbb_eur.tsv', index=False, header=False, sep='\t')

high_missingness_low_maf_variants = concat([high_missingness_variants, low_maf_variants]).drop_duplicates(ignore_index=True)

high_missingness_low_maf_variants.to_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/high_missingness_low_maf_variants_panukbb_eur.tsv', index=False, header=False, sep='\t')

# 7. Calculate HWE stats

# 8. Remove variants based on HWE stats
disease_fields = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')['disease_field'].to_list()
num_vars_excluded_by_hwe_list = []
bar = ProgBar(len(disease_fields), stream=sys.stdout, title='Generating variants to exclude')
for disease in disease_fields:
    autosomes_hwe = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step7_hwe/' + disease + '_panukbb_eur.hardy', delim_whitespace=True)
    autosomes_hwe_fail = autosomes_hwe[autosomes_hwe['MIDP'] < hwe_autosomes_threshold]

    cx_hwe = read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step7_hwe/' + disease + '_panukbb_eur.hardy.x', delim_whitespace=True)
    cx_hwe_fail = cx_hwe[cx_hwe['MIDP'] < hwe_autosomes_threshold]

    hwe_fail_variants = concat([autosomes_hwe_fail[['ID']], cx_hwe_fail[['ID']]]).drop_duplicates(ignore_index=True)
    num_vars_excluded_by_hwe_list.append(len(hwe_fail_variants))

    concat([high_missingness_low_maf_variants, hwe_fail_variants]).to_csv(
        '~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/variants_to_exclude/' + disease + '_autosomes-cX-vars_to_exclude.tsv',
        header=False, index=False, sep='\t')

    bar.update()
num_vars_excluded_by_hwe_per_disease = DataFrame({'disease_field': disease_fields,
                                                  'num_vars_excluded_by_hwe': num_vars_excluded_by_hwe_list})
num_vars_excluded_by_hwe_per_disease.to_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step8_num_vars_excluded_by_hwe.csv', index=False)

""" Pipeline before 2024-02-10 (Genomic QC 1)"""

# # minor allele frequency
# autosomes_freq = read_csv('ukbiobank/gwas/snp_qc/stats/autosomes_panukbb_eur.frq', delim_whitespace=True)
# autosomes_freq_fail = autosomes_freq[autosomes_freq['MAF'] < maf_threshold]
#
# cx_males_freq = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_male.frq', delim_whitespace=True)
# cx_females_freq = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_female.frq', delim_whitespace=True)
# cx_freq_fail = concat([cx_males_freq[cx_males_freq['MAF'] < maf_threshold],
#                        cx_females_freq[cx_females_freq['MAF'] < maf_threshold]])
#
# # missingness: variants
# autosomes_lmiss = read_csv('ukbiobank/gwas/snp_qc/stats/autosomes_panukbb_eur.lmiss', delim_whitespace=True)
# autosomes_lmiss['F_MISS_CALCED'] = (autosomes_lmiss['N_MISS'] / autosomes_lmiss['N_GENO'])
# autosomes_lmiss_fail = autosomes_lmiss[autosomes_lmiss['F_MISS_CALCED'] > lmiss_autosomes_threshold]
#
# cx_males_lmiss = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_male.lmiss', delim_whitespace=True)
# cx_males_lmiss['F_MISS_CALCED'] = (cx_males_lmiss['N_MISS'] / cx_males_lmiss['N_GENO'])
# cx_females_lmiss = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_female.lmiss', delim_whitespace=True)
# cx_females_lmiss['F_MISS_CALCED'] = (cx_females_lmiss['N_MISS'] / cx_females_lmiss['N_GENO'])
# cx_mf_diff_lmiss = DataFrame({'SNP': cx_females_lmiss['SNP'],
#                              'male_female_diff': abs(cx_females_lmiss['F_MISS_CALCED'] - cx_males_lmiss['F_MISS_CALCED'])})
# cx_lmiss_fail = concat([cx_females_lmiss[cx_females_lmiss['F_MISS_CALCED'] > lmiss_cx_threshold],
#                         cx_males_lmiss[cx_males_lmiss['F_MISS_CALCED'] > lmiss_cx_threshold],
#                         cx_mf_diff_lmiss[cx_mf_diff_lmiss['male_female_diff'] > lmiss_cx_threshold]])
#
# # missingness: samples
# autosomes_imiss = read_csv('ukbiobank/gwas/snp_qc/stats/autosomes_panukbb_eur.imiss', delim_whitespace=True)
# autosomes_imiss['F_MISS_CALCED'] = (autosomes_imiss['N_MISS'] / autosomes_imiss['N_GENO'])
# autosomes_imiss_fail = autosomes_imiss[autosomes_imiss['F_MISS_CALCED'] > imiss_threshold]
#
# cx_males_imiss = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_male.imiss', delim_whitespace=True)
# cx_males_imiss['F_MISS_CALCED'] = (cx_males_imiss['N_MISS'] / cx_males_imiss['N_GENO'])
# cx_females_imiss = read_csv('ukbiobank/gwas/snp_qc/stats/cX_panukbb_eur_female.imiss', delim_whitespace=True)
# cx_females_imiss['F_MISS_CALCED'] = (cx_females_imiss['N_MISS'] / cx_females_imiss['N_GENO'])
# cx_imiss_fail = concat([cx_males_imiss[cx_males_imiss['F_MISS_CALCED'] > imiss_threshold],
#                         cx_females_imiss[cx_females_imiss['F_MISS_CALCED'] > imiss_threshold]])
#
# # concat failed variants
# autosome_variants_to_remove = concat([autosomes_freq_fail[['SNP']],
#                                       autosomes_lmiss_fail[['SNP']]]).drop_duplicates(ignore_index=True)
# cx_variants_to_remove = concat([cx_freq_fail[['SNP']],
#                                 cx_lmiss_fail[['SNP']]]).drop_duplicates(ignore_index=True)
#
# # filter out, per disease, variants that fail HWE QC, and append to the variants that failed the MAF/missingness QC,
# # to generate disease-specific lists of variants to be excluded
# disease_fields = read_csv('ukbiobank/ard_identification/by_age_of_onset_stats/qu10_50_diseases_icd10-info_num-diagnosed-info.csv')['disease_field'].to_list()
# bar = ProgBar(len(disease_fields), stream=sys.stdout, title='Generating variants to exclude')
# for disease in disease_fields:
#
#     autosomes_hwe = read_csv('ukbiobank/gwas/snp_qc/stats/hwe/' + disease + '_hwe.hardy', delim_whitespace=True)
#     autosomes_hwe_fail = autosomes_hwe[autosomes_hwe['MIDP'] < hwe_autosomes_threshold]
#     autosome_variants_to_remove_disease = concat([autosome_variants_to_remove,
#                                                   autosomes_hwe_fail[['ID']]]).drop_duplicates(ignore_index=True)
#     autosome_variants_to_remove_disease.to_csv(
#         'ukbiobank/gwas/snp_qc/variants_to_exclude/' + disease + '_autosome-vars_to_exclude.tsv', header=False,
#         index=False, sep='\t')
#
#     cx_hwe = read_csv('ukbiobank/gwas/snp_qc/stats/hwe/' + disease + '_hwe.hardy.x', delim_whitespace=True)
#     cx_hwe_fail = cx_hwe[cx_hwe['MIDP'] < hwe_autosomes_threshold]
#     cx_variants_to_remove_disease = concat([cx_variants_to_remove,
#                                             cx_hwe_fail[['ID']]]).drop_duplicates(ignore_index=True)
#     cx_variants_to_remove_disease.to_csv('ukbiobank/gwas/snp_qc/variants_to_exclude/' + disease + '_cX-vars_to_exclude.tsv',
#                                          header=False, index=False, sep='\t')
#
#     variants_to_remove = concat([autosome_variants_to_remove_disease, cx_variants_to_remove_disease])
#     variants_to_remove.to_csv('ukbiobank/gwas/snp_qc/variants_to_exclude/' + disease + '_autosomes-cX-vars_to_exclude.tsv',
#                               header=False, index=False, sep='\t')
#
#     bar.update()
