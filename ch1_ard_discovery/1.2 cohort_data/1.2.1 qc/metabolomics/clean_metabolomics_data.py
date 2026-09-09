"""
Script to clean the UKB NMR metabolomics data downloaded off DNAnexus. Following tasks performed:
- Determine sex and Pan-UK Biobank ancestry group assignments for all participants
- Isolate EIDs for Pan-UKBB EUR
- Isolate all metabolomics data (covars, QC flags, data) for Pan-UKBB EUR
- Sequentially remove metabolites and participants with missingness > 5%
- log1p-transform data
- Separate whole data into constituent parts (covars, QC flags, data, eids) and export
"""
from pandas import read_csv, DataFrame, merge
from pyprind import ProgBar
from datetime import datetime

import numpy as np
import sys

metabolomics = read_csv('~/data/external/metabolomics/metabolomics_data_covars_qcflags_off_dnanexus.csv')

data_fields = ['p' + str(x) + '_i0' for x in list(range(23400, 23649))]

""" Sex and ancestry group assignment to the 'everyone' from the metabolomics cohort """
pan_ukbb_ancestry = read_csv('~/data/external/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_eid98294.tsv')
pan_ukbb_ancestry_per_person = pan_ukbb_ancestry[['eid98294', 'pop']]
pan_ukbb_ancestry_per_person.columns = ['eid', 'pan_ukbb_ancestry_group']

eid_sex = metabolomics[['eid', 'p31']]
eid_sex.columns = ['eid', 'sex']

eid_sex_ancestry = merge(eid_sex, pan_ukbb_ancestry_per_person, on='eid')
# eid_sex_ancestry.to_csv('ukbiobank/by_sex_ethnicity/sex-ancestry-group_per_metabolome_person.csv', index=False)

""" Isolate Pan-UKBB EUR """
pan_ukbb_eur_eids = eid_sex_ancestry[eid_sex_ancestry.pan_ukbb_ancestry_group == 'EUR'][['eid']]
metabolomics_eur = metabolomics[metabolomics.eid.isin(list(pan_ukbb_eur_eids['eid']))].copy()

""" Isolate those also present in the genomics cohort """
pan_ukbb_eur_eids_genomics = list(read_csv('~/data/internal/cohort_eids/genomics_pan-ukbb-eur_eids.tsv', delim_whitespace=True, header=None)[0])
high_missingness_eids_genomics = list(read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
                                      delim_whitespace=True, header=None)[0])
metabolomics_eur_genetic = metabolomics_eur[(metabolomics_eur['eid'].isin(pan_ukbb_eur_eids_genomics)) & (~metabolomics_eur['eid'].isin(high_missingness_eids_genomics))]

""" Assess missingness """
metabolomics_data_raw = metabolomics_eur_genetic.loc[:, ['eid'] + data_fields]
metabolite_missingness_pre = metabolomics_data_raw.iloc[:, 1:].isna().sum() / len(metabolomics_data_raw)
sample_missingness_pre = metabolomics_data_raw.iloc[:, 1:].isna().sum(axis=1) / len(data_fields)

""" Filter metabolites by missingness > 5% """
high_missingness_metabolites = list(metabolite_missingness_pre[metabolite_missingness_pre > 0.05].index)
metabolomics_data_no_high_missingness_metabolites = metabolomics_data_raw.drop(high_missingness_metabolites, axis=1)

""" Filter samples by missingness > 5% """
sample_missingness_no_high_missingness_metabolites = metabolomics_data_no_high_missingness_metabolites.iloc[
                                                     :, 1:].isna().sum(axis=1) / len(data_fields)
high_missingness_samples = list(sample_missingness_no_high_missingness_metabolites[
                                    sample_missingness_no_high_missingness_metabolites > 0.05].index)
metabolomics_data_post_missingness = metabolomics_data_no_high_missingness_metabolites.drop(high_missingness_samples)

""" Re-assess missingness """
metabolite_missingness_post = metabolomics_data_post_missingness.iloc[:, 1:].isna().sum() / len(metabolomics_data_post_missingness)
sample_missingness_post = metabolomics_data_post_missingness.iloc[:, 1:].isna().sum(axis=1) / (len(data_fields) - len(high_missingness_metabolites))

""" log1p-transforming data """
metabolomics_data_post_missingness_log1p = metabolomics_data_post_missingness.copy()
metabolomics_data_post_missingness_log1p.iloc[:, 1:] = np.log1p(metabolomics_data_post_missingness.iloc[:, 1:])

""" Export eid + QC data """
first_qc_field_col_index = metabolomics_eur_genetic.columns.get_loc(data_fields[-1]) + 1
qc_fields_indices = list(range(first_qc_field_col_index, metabolomics_eur_genetic.shape[1]))
metabolomics_qc_flags = metabolomics_eur_genetic.iloc[:, [0] + qc_fields_indices].drop(high_missingness_samples)
for metabolite in high_missingness_metabolites:
    metabolite_qc_field = 'p' + str(int(metabolite[1:6]) + 300) + '_i0'
    if metabolite_qc_field in metabolomics_qc_flags.columns:
        metabolomics_qc_flags.drop(metabolite_qc_field, axis=1, inplace=True)
metabolomics_qc_flags.to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_qc_flags.csv', index=False)

""" Export eid + covariates """
metabolomics_covars = metabolomics_eur_genetic.loc[:, :data_fields[0]].drop(
    [data_fields[0]] + high_missingness_metabolites, axis=1).drop(high_missingness_samples)

# calculate storage time
times_between_sample_taken_and_measurement = []
for i in range(len(metabolomics_covars)):
    date_sample_taken = datetime.strptime(metabolomics_covars.loc[i, 'p53_i0'], '%Y-%m-%d')
    date_sample_measured = datetime.fromisoformat(metabolomics_covars.loc[i, 'p23659_i0'])
    timezone = date_sample_measured.tzinfo
    date_sample_taken = date_sample_taken.replace(tzinfo=timezone)
    times_between_sample_taken_and_measurement.append((date_sample_measured - date_sample_taken).days)
metabolomics_covars['storage_time'] = times_between_sample_taken_and_measurement

# calculate time between sample prep and measurement
times_between_sample_prep_and_measurement = []
for i in range(len(metabolomics_covars)):
    times_between_sample_prep_and_measurement.append((datetime.fromisoformat(metabolomics_covars.loc[i, 'p23658_i0']) - datetime.fromisoformat(metabolomics_covars.loc[i, 'p23659_i0'])).total_seconds())
metabolomics_covars['prepped_for_time'] = times_between_sample_prep_and_measurement

metabolomics_covars['plate_row'] = metabolomics_covars['p23660_i0'].apply(lambda x: x[0])
metabolomics_covars['plate_column'] = metabolomics_covars['p23660_i0'].apply(lambda x: x[1:])

metabolomics_covars['sample_prepared_date'] = metabolomics_covars['p23659_i0'].apply(lambda x: datetime.fromisoformat(x).strftime('%Y-%m-%d'))
metabolomics_covars['sample_measured_date'] = metabolomics_covars['p23658_i0'].apply(lambda x: datetime.fromisoformat(x).strftime('%Y-%m-%d'))

metabolomics_covars.to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_covars.csv', index=False)

""" Export eid + NMR data """
# metabolomics_data_post_missingness.to_csv('ukbiobank/metabolomics/metabolomics_pan_ukbb_eur_data.csv', index=False)

""" Export eid + log1p-transformed NMR data """
metabolomics_data_post_missingness_log1p.to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_log1p_data.csv', index=False)

""" Export eid """
metabolomics_data_post_missingness[['eid']].to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_eids.csv', index=False)

""" Export eid + lifestyle factors """
metabolomics_lifestyle = read_csv('~/data/external/metabolomics/metabolomics_lifestyles_off_dnanexus.csv')
metabolomics_lifestyle_eur = metabolomics_lifestyle[metabolomics_lifestyle['eid'].isin(list(metabolomics_data_post_missingness.eid))].copy()
# lifestyle_nas = metabolomics_lifestyle_eur.isna().sum()
# indiv_nas = metabolomics_lifestyle_eur.isna().sum(axis=1)
# indiv_nas_many = indiv_nas[indiv_nas > 5].sort_values()
# codify certain values and impute missing values
# alcohol - categorical
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p1558_i0.isna(), 'p1558_i0'] = 'Prefer not to answer'
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p1558_i0 == 'Prefer not to answer', 'p1558_i0'] = 'Unanswered'
# smoking status - categorical
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p20116_i0.isna(), 'p20116_i0'] = 'Prefer not to answer'
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p20116_i0 == 'Prefer not to answer', 'p1558_i0'] = 'Unanswered'
# % body fat - continouous
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p23099_i0.isna(), 'p23099_i0'] = np.nanmedian(metabolomics_lifestyle_eur['p23099_i0'])
# hours of sleep - continuous
metabolomics_lifestyle_eur.loc[~metabolomics_lifestyle_eur['p1160_i0'].isin([str(x) for x in range(24)]), 'p1160_i0'] = round(np.nanmedian(metabolomics_lifestyle_eur[metabolomics_lifestyle_eur['p1160_i0'].isin([str(x) for x in range(24)])]['p1160_i0'].astype(int)))
metabolomics_lifestyle_eur.loc[:, 'p1160_i0'] = metabolomics_lifestyle_eur.p1160_i0.astype(int)
# Townsend deprivation index - continuous
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p22189.isna(), 'p22189'] = np.nanmedian(metabolomics_lifestyle_eur['p22189'])
# education - categorical
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur['p6138_i0'].isna(), 'p6138_i0'] = 'no_higher_education'
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur['p6138_i0'].str.contains('College or University degree'), 'p6138_i0'] = 'higher_education'
metabolomics_lifestyle_eur.loc[~metabolomics_lifestyle_eur['p6138_i0'].str.contains('higher_education'), 'p6138_i0'] = 'no_higher_education'
# physical activity - continuous
metabolomics_lifestyle_eur.loc[~metabolomics_lifestyle_eur['p884_i0'].isin([str(x) for x in range(8)]), 'p884_i0'] = round(np.nanmedian(metabolomics_lifestyle_eur[metabolomics_lifestyle_eur['p884_i0'].isin([str(x) for x in range(8)])]['p884_i0'].astype(int)))
metabolomics_lifestyle_eur.loc[:, 'p884_i0'] = metabolomics_lifestyle_eur.p884_i0.astype(int)
# grip strength - continuous
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p46_i0.isna(), 'p46_i0'] = np.nanmedian(metabolomics_lifestyle_eur.p46_i0)
metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur.p47_i0.isna(), 'p47_i0'] = np.nanmedian(metabolomics_lifestyle_eur.p47_i0)
metabolomics_lifestyle_eur['grip_strength'] = (metabolomics_lifestyle_eur.p46_i0 + metabolomics_lifestyle_eur.p47_i0) / 2
# various food items - continuous
for field in ['p1289_i0', 'p1299_i0', 'p1309_i0', 'p1319_i0', 'p1438_i0', 'p1458_i0', 'p1528_i0']:
    metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur[field].isna(), field] = 'Do not know'
    metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur[field] == 'Less than one', field] = '0.5'
    metabolomics_lifestyle_eur.loc[
        metabolomics_lifestyle_eur[field].isin(['Prefer not to answer', 'Do not know']), field] = round(
        np.nanmedian(metabolomics_lifestyle_eur[metabolomics_lifestyle_eur[field].isin([str(x) for x in range(101)])][
                         field].astype(float)))
    metabolomics_lifestyle_eur.loc[:, field] = metabolomics_lifestyle_eur[field].astype(float)
# various food items - categorical
for field in ['p' + str(1329 + 10*x) + '_i0' for x in range(7)]:
    metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur[field].isna(), field] = 'Do not know'
    metabolomics_lifestyle_eur.loc[metabolomics_lifestyle_eur[field].isin(['Prefer not to answer', 'Do not know']), field] = 'Unanswered'
# export
metabolomics_lifestyle_eur.drop(['p46_i0', 'p47_i0', 'p48_i0', 'p49_i0'], axis=1).to_csv('~/data/internal/metabolomics_pan_ukbb_eur_lifestyles.csv', index=False)

""" Note data fields which are ratios """
metabolomics_fields = read_csv('~/data/external/nmr_nightingale_metabolites_grouped_listing.txt', sep='\t')
metabolomics_fields.fillna('', inplace=True)
metabolomics_ratio_fields = metabolomics_fields[(metabolomics_fields['title'].str.contains(' ratio ')) |
                                                (metabolomics_fields['title'].str.contains('Ratio '))].copy()
metabolomics_ratio_fields.loc[:, 'field'] = metabolomics_ratio_fields['field_id'].apply(lambda x: 'p' + str(x) + '_i0')
metabolomics_ratio_fields.to_csv('~/data/internal/metabolomics/metabolomics_ratio_data_fields.csv', index=False)

""" Filter phenotype data to new eid set """
metabolomics_phenotypes = read_csv('~/data/external/metabolomics/metabolomics_phenotype_file_off_dnanexus.csv', low_memory=True)
metabolomics_phenotypes = metabolomics_phenotypes[metabolomics_phenotypes['eid'].isin(list(metabolomics_covars['eid']))]
metabolomics_phenotypes.to_csv('~/data/internal/metabolomics/metabolomics/metabolomics_pan_ukbb_eur_phenotypes.csv', index=False)
