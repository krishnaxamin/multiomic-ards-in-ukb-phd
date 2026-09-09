"""
Extract num_diagnosed information for common-unisex diseases from different cohorts.
"""
import pandas as pd

""" Genomics """

phenotypes = pd.read_csv('~/data/external/genomics/array-genotyping_phenotype_file.csv')

pan_ukbb_eur_eids = list(
    pd.read_csv('~/data/internal/cohort_eids/genomics_pan-ukbb-eur_eids.tsv', delim_whitespace=True, header=None)[0])
high_missingness_eids = list(pd.read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
                                         delim_whitespace=True, header=None)[0])
filtered_phenotypes = phenotypes[
    (phenotypes['eid'].isin(pan_ukbb_eur_eids)) & (~phenotypes['eid'].isin(high_missingness_eids))]
common_unisex_fields_df = pd.read_csv('~/data/internal/genomics/genomics_qc2_pan_ukbb_eur_common_unisex_disease_date_fields.csv')
common_unisex_fields = list(common_unisex_fields_df[common_unisex_fields_df['everyone'] == 1]['disease_field'])

num_in_cohort = len(filtered_phenotypes)

field_dfs = []
for field in common_unisex_fields:
    num_diagnosed = filtered_phenotypes[field].notnull().sum()
    case_control_ratio = (num_in_cohort - num_diagnosed) / num_diagnosed
    field_df = pd.DataFrame({'disease_field': [field],
                             'num_cases': [num_diagnosed],
                             'case_control_ratio': [case_control_ratio]})
    field_dfs.append(field_df)

total_df = pd.concat(field_dfs)
total_df.to_csv('~/data/internal/genomics/genomics_qc2_pan_ukbb_eur_common-unisex-disease_num-diagnosed_case-control.csv', index=False)

""" Metabolomics """

metabolomics_phenotypes = pd.read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_phenotypes.csv', low_memory=False)
common_unisex_fields_df = pd.read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_common_unisex_disease_date_fields.csv')
common_unisex_fields = list(common_unisex_fields_df[common_unisex_fields_df['everyone'] == 1]['disease_field'])

num_in_cohort = len(metabolomics_phenotypes)

field_dfs = []
for field in common_unisex_fields:
    num_diagnosed = metabolomics_phenotypes[field].notnull().sum()
    case_control_ratio = (num_in_cohort - num_diagnosed) / num_diagnosed
    field_df = pd.DataFrame({'disease_field': [field],
                             'num_cases': [num_diagnosed],
                             'case_control_ratio': [case_control_ratio]})
    field_dfs.append(field_df)

total_df = pd.concat(field_dfs)
total_df.to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_common-unisex-disease_num-diagnosed_case-control.csv', index=False)

""" Proteomics """

proteomics_phenotypes = pd.read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_phenotypes.csv', low_memory=False)
common_unisex_fields_df = pd.read_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_common_unisex_disease_date_fields.csv')
common_unisex_fields = list(common_unisex_fields_df[common_unisex_fields_df['everyone'] == 1]['disease_field'])

num_in_cohort = len(proteomics_phenotypes)

field_dfs = []
for field in common_unisex_fields:
    num_diagnosed = proteomics_phenotypes[field].notnull().sum()
    case_control_ratio = (num_in_cohort - num_diagnosed) / num_diagnosed
    field_df = pd.DataFrame({'disease_field': [field],
                             'num_cases': [num_diagnosed],
                             'case_control_ratio': [case_control_ratio]})
    field_dfs.append(field_df)

total_df = pd.concat(field_dfs)
total_df.to_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_common-unisex-disease_num-diagnosed_case-control.csv',
                index=False)
