""" LCV results' relationship to ICD-10 chapters. """
import scipy
import pandas as pd

lcv = pd.read_csv('~/data/internal/genomics/causality/lcv_results_processed.csv')

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

lcv = (lcv
       .merge(disease_info[['disease_field', 'code_chapter']]
              .rename(columns={'disease_field': 'disease1'}))
       .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
       .merge(disease_info[['disease_field', 'code_chapter']]
              .rename(columns={'disease_field': 'disease2'}))
       .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))

lcv_flagged = lcv[(lcv['gcp0_fdr'] < 0.05) & (abs(lcv['gcp_est']) > 0.6)].copy()
lcv_flagged_slim = lcv_flagged[['disease1', 'disease2', 'gcp_est']].copy()

""" Relationship to ICD-10 chapters """
lcv_flagged_slim_chapters = (lcv_flagged_slim
                  .merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease1', 'icd10_chapter': 'disease1_chapter'}))
                  .merge(disease_info[['code_chapter', 'icd10_chapter']].rename(columns={'code_chapter': 'disease2', 'icd10_chapter': 'disease2_chapter'}))
                  .drop_duplicates())
lcv_flagged_slim_chapters['chapter_match'] = lcv_flagged_slim_chapters['disease1_chapter'] == lcv_flagged_slim_chapters['disease2_chapter']
X = lcv_flagged_slim_chapters['chapter_match'].sum()
Y = len(lcv_flagged_slim_chapters)
scipy.stats.fisher_exact([[X, 293 - X], [Y - X, 1985 - Y + X]], alternative='greater')

# GCP in same-chapter pairs vs diff-chapter pairs
scipy.stats.mannwhitneyu(abs(lcv_flagged_slim_chapters[lcv_flagged_slim_chapters['chapter_match'] == 1].gcp_est),
                         abs(lcv_flagged_slim_chapters[lcv_flagged_slim_chapters['chapter_match'] == 0].gcp_est))

# GCP in concordant pairs vs discordant pairs
scipy.stats.mannwhitneyu(abs(lcv_flagged_slim_chapters[lcv_flagged_slim_chapters['lcv_mw_concordant'] == 1].gcp_est),
                         abs(lcv_flagged_slim_chapters[lcv_flagged_slim_chapters['lcv_mw_concordant'] == 0].gcp_est))
