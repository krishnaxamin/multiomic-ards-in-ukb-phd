""" Compare LCV results with directional disease relationships based on relative age-of-onset. """
from pandas import read_csv, merge, concat, DataFrame
from scipy.stats import mannwhitneyu
from itertools import chain

disease_info = read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

lcv = read_csv('~/data/internal/genomics/causality/lcv_results_processed.csv')

# validate on age-of-onset distributions
#   disease1 -> disease2 should have disease2 having later age-of-onset than disease1
lcv_flagged = lcv[(lcv['gcp0_fdr'] < 0.05) & (abs(lcv['gcp_est']) > 0.6)].copy()
lcv_flagged_slim = lcv_flagged[['disease1', 'disease2', 'gcp_est']].reset_index(
    drop=True).copy()

mw_results = []
for df_row in lcv_flagged_slim.to_dict(orient='records'):
    disease1_time_series = read_csv(
        '~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/' + df_row[
            'disease1'] + '_everyone_time_series.csv')
    disease1_age_set = list(chain.from_iterable(
        [[dic['age']] * int(dic['num_diagnosed']) for dic in disease1_time_series.to_dict(orient='records')]))
    disease2_time_series = read_csv(
        '~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/' + df_row[
            'disease2'] + '_everyone_time_series.csv')
    disease2_age_set = list(chain.from_iterable(
        [[dic['age']] * int(dic['num_diagnosed']) for dic in disease2_time_series.to_dict(orient='records')]))
    mw_results.append({
        'disease1': df_row['disease1'],
        'disease2': df_row['disease2'],
        'mw_pval_two_sided': mannwhitneyu(disease1_age_set, disease2_age_set).pvalue,
        'mw_pval_greater': mannwhitneyu(disease1_age_set, disease2_age_set,
                                        alternative='greater').pvalue,
        'mw_pval_less': mannwhitneyu(disease1_age_set, disease2_age_set,
                                     alternative='less').pvalue
    })
mw_results_df = DataFrame(mw_results)
lcv_mw_df = lcv_flagged_slim.merge(mw_results_df)
lcv_mw_df['lcv_mw_concordant'] = 0
lcv_mw_df.loc[(((lcv_mw_df['gcp_est'] > 0) & (lcv_mw_df['mw_pval_less'] < 0.05)) |
               ((lcv_mw_df['gcp_est'] < 0) & (lcv_mw_df['mw_pval_greater'] < 0.05))), 'lcv_mw_concordant'] = 1
lcv_mw_df_codes = (lcv_mw_df
                   .merge(disease_info[['disease_field', 'code_chapter']]
                          .rename(columns={'disease_field': 'disease1'}))
                   .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                   .merge(disease_info[['disease_field', 'code_chapter']]
                          .rename(columns={'disease_field': 'disease2'}))
                   .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))
lcv_mw_df_codes.to_csv('~/data/internal/genomics/causality/lcv_concordance_with_age_of_onset_fdr0.05_gcp0.6.csv', index=False)

""" Relationship to ICD-10 chapters """
# concordant set enriched in same-chapter pairs?
import scipy
import pandas as pd
lcv_mw_df_codes = pd.read_csv('~/data/internal/genomics/causality/lcv_concordance_with_age_of_onset_fdr0.05_gcp0.6.csv')
lcv_concordant = lcv_mw_df_codes[lcv_mw_df_codes['lcv_mw_concordant'] == 1][['disease1', 'disease2']].copy()
lcv_concordant = (lcv_concordant
                  .merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease1', 'icd10_chapter': 'disease1_chapter'}))
                  .merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease2', 'icd10_chapter': 'disease2_chapter'}))
                  .drop_duplicates())
lcv_concordant['chapter_match'] = lcv_concordant['disease1_chapter'] == lcv_concordant['disease2_chapter']
X = lcv_concordant['chapter_match'].sum()
Y = len(lcv_concordant)
scipy.stats.fisher_exact([[X, 293 - X], [Y - X, 1985 - Y + X]], alternative='greater')

# GCP in same-chapter pairs vs diff-chapter pairs
scipy.stats.mannwhitneyu(abs(lcv_concordant[lcv_concordant['chapter_match'] == 1].gcp_est),
                         abs(lcv_concordant[lcv_concordant['chapter_match'] == 0].gcp_est))

# GCP in concordant pairs vs discordant pairs
scipy.stats.mannwhitneyu(abs(lcv_mw_df_codes[lcv_mw_df_codes['lcv_mw_concordant'] == 1].gcp_est),
                         abs(lcv_mw_df_codes[lcv_mw_df_codes['lcv_mw_concordant'] == 0].gcp_est))
