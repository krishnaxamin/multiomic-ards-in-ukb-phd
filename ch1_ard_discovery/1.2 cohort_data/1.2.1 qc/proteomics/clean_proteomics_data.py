"""
Script to clean the UKB proteomics data downloaded off DNAnexus. Following tasks performed:
- Check participants do not include those selected by UKB-PPP consortium
- Determine sex and Pan-UK Biobank ancestry group assignments for all participants
- Isolate EIDs for Pan-UKBB EUR
- Isolate all proteomics for Pan-UKBB EUR
- Remove data from batches 0 and 7
- Sequentially remove proteins and participants with missingness > 25%
- Merge in additional covariate data from Olink data files, setting up participant- and protein-based covariate files
- Export data in constituent parts (covars, data, eids)
"""
from pandas import read_csv, merge, DataFrame, concat
from pyprind import ProgBar

import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt

""" Load data """
npx = read_csv('~/data/external/proteomics/proteomics_npx_data_off_dnanexus.csv')
non_npx = read_csv('~/data/external/proteomics/proteomics_non_npx_data_off_dnanexus.csv')

olink_lod = read_csv('~/data/external/proteomics/olink_limit_of_detection.dat', delim_whitespace=True)
olink_assay_warning = read_csv('~/data/external/proteomics/olink_assay_warning.dat', delim_whitespace=True)  # all PASS
olink_assay = read_csv('~/data/external/proteomics/olink_assay.dat', sep='\t')
olink_batch = read_csv('~/data/external/proteomics/olink_batch_number.dat', delim_whitespace=True)
olink_panel_lot = read_csv('~/data/external/proteomics/olink_panel_lot_number.dat', sep='\t')
olink_assay_version = read_csv('~/data/external/proteomics/olink_assay_version.dat',
                               delim_whitespace=True)  # all version 1
olink_processing_start_date = read_csv('~/data/external/proteomics/olink_processing_start_date.dat', sep='\t')

""" Check participants are not those selected by UKB-PPP consortium """
non_npx.p30903_i0.notna().sum()

""" Sex and ancestry group assignment to the 'everyone' from the proteomics cohort """
pan_ukbb_ancestry = read_csv('~/data/external/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_eid98294.tsv')
pan_ukbb_ancestry_per_person = pan_ukbb_ancestry[['eid98294', 'pop']]
pan_ukbb_ancestry_per_person.columns = ['eid', 'pan_ukbb_ancestry_group']

eid_sex = non_npx[['eid', 'p31']]
eid_sex.columns = ['eid', 'sex']

eid_sex_ancestry = merge(eid_sex, pan_ukbb_ancestry_per_person, on='eid')
# eid_sex_ancestry.to_csv('ukbiobank/by_sex_ethnicity/sex-ancestry-group_per_proteome_person.csv', index=False)

""" Isolate Pan-UKBB EUR """
pan_ukbb_eur_eids = eid_sex_ancestry[eid_sex_ancestry.pan_ukbb_ancestry_group == 'EUR'][['eid']]
npx_eur = npx[npx.eid.isin(list(pan_ukbb_eur_eids['eid']))].copy()

""" Isolate those also preset in the genomics cohort """
pan_ukbb_eur_eids_genomics = list(
    read_csv('~/data/internal/cohort_eids/genomics_pan-ukbb-eur_eids.tsv', delim_whitespace=True, header=None)[0])
high_missingness_eids_genomics = list(
    read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
             delim_whitespace=True, header=None)[0])
npx_eur_genetic = npx_eur[
    (npx_eur['eid'].isin(pan_ukbb_eur_eids_genomics)) & (~npx_eur['eid'].isin(high_missingness_eids_genomics))]

""" Add plate info to NPX data """
npx_plate = merge(npx_eur_genetic, non_npx[['eid', 'p30901_i0']])
npx_plate.rename(columns={'p30901_i0': 'PlateID'}, inplace=True)

""" Remove data from batches 0 and 7 """
olink_plates_in_batches_1_to_6 = list(olink_batch[(olink_batch['Batch'] != 0) & (olink_batch['Batch'] != 7)]['PlateID'])
npx_plate_batch_filtered = npx_plate[npx_plate['PlateID'].isin(olink_plates_in_batches_1_to_6)].copy().reset_index(
    drop=True)

# total number of data points after filtering for batches
num_npx_data_points_batch_filtered = npx_plate_batch_filtered.iloc[:, 1:-1].notna().sum().sum()

""" Calculate and filter protein missingness """
protein_missingness = npx_plate_batch_filtered.iloc[:, 1:-1].isna().sum() / len(npx_plate_batch_filtered)
# (protein_missingness > 0.25).sum()
high_missingness_proteins = list(protein_missingness[protein_missingness > 0.25].index)
DataFrame({'protein_field': high_missingness_proteins}).to_csv(
    '~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/proteomics/high_missingness_protein_fields.txt',
    index=False)
npx_plate_batch_filtered_no_high_missingness_proteins = npx_plate_batch_filtered.drop(high_missingness_proteins, axis=1)

""" Calculate and filter on sample missingness """
all_proteins = list(npx_eur_genetic.columns[1:])
sample_missingness = npx_plate_batch_filtered_no_high_missingness_proteins.iloc[:, 1:].isna().sum(axis=1) / (
            len(npx_plate_batch_filtered_no_high_missingness_proteins.columns) - 2)
high_missingness_samples = list(sample_missingness[sample_missingness > 0.25].index)
npx_plate_batch_filtered_no_high_missingness_proteins_samples = npx_plate_batch_filtered_no_high_missingness_proteins.drop(
    high_missingness_samples)

""" Merge in more covariate data from Olink data files """
covars = merge(non_npx, npx_plate_batch_filtered_no_high_missingness_proteins_samples[['eid']],
               on='eid')  # filter to correct EIDs
covars = merge(covars, olink_batch.rename(columns={'PlateID': 'p30901_i0'}), on='p30901_i0')  # add batch
# separate well position into well column and well row, as done in Ritchie et al, 2023
covars['plate_row'] = covars['p30902_i0'].apply(lambda x: x[0])
covars['plate_column'] = covars['p30902_i0'].apply(lambda x: x[1:])

protein_based_covars = merge(olink_assay, olink_panel_lot.drop('Batch', axis=1), on='Panel').drop_duplicates()
protein_based_covars = protein_based_covars[
    ~protein_based_covars['Assay'].isin([x.upper() for x in high_missingness_proteins])].copy()
protein_based_covars = merge(protein_based_covars, olink_processing_start_date, on='Panel')
protein_based_covars = protein_based_covars[protein_based_covars['PlateID'].isin(list(set(covars['p30901_i0'])))].copy()

""" Export eid + covariates """
# covars = merge(non_npx, npx_plate_batch_filtered_no_high_missingness_proteins_samples[['eid']], on='eid')
covars.to_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv', index=False)
protein_based_covars.to_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv', index=False)

""" Export eid + proteomics data """
npx_plate_batch_filtered_no_high_missingness_proteins_samples.drop('PlateID', axis=1).sort_values(by='eid').to_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_data.csv', index=False)

""" Export eid """
npx_plate_batch_filtered_no_high_missingness_proteins_samples[['eid']].to_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_eids.csv', index=False)

""" Export eid + lifestyle factors """
proteomics_lifestyle = read_csv('ukbiobank/proteomics/proteomics_lifestyles_off_dnanexus.csv')
proteomics_lifestyle_eur = proteomics_lifestyle[
    proteomics_lifestyle['eid'].isin(list(npx_plate_batch_filtered_no_high_missingness_proteins_samples.eid))].copy()
lifestyle_nas = proteomics_lifestyle_eur.isna().sum()
# indiv_nas = proteomics_lifestyle_eur.isna().sum(axis=1)
# indiv_nas_many = indiv_nas[indiv_nas > 5].sort_values()
# codify certain values and impute missing values
# alcohol - categorical
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p1558_i0.isna(), 'p1558_i0'] = 'Prefer not to answer'
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p1558_i0 == 'Prefer not to answer', 'p1558_i0'] = 'Unanswered'
# smoking status - categorical
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p20116_i0.isna(), 'p20116_i0'] = 'Prefer not to answer'
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p20116_i0 == 'Prefer not to answer', 'p1558_i0'] = 'Unanswered'
# % body fat - continouous
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p23099_i0.isna(), 'p23099_i0'] = np.nanmedian(
    proteomics_lifestyle_eur['p23099_i0'])
# hours of sleep - continuous
proteomics_lifestyle_eur.loc[
    ~proteomics_lifestyle_eur['p1160_i0'].isin([str(x) for x in range(24)]), 'p1160_i0'] = round(np.nanmedian(
    proteomics_lifestyle_eur[proteomics_lifestyle_eur['p1160_i0'].isin([str(x) for x in range(24)])]['p1160_i0'].astype(
        int)))
proteomics_lifestyle_eur.loc[:, 'p1160_i0'] = proteomics_lifestyle_eur.p1160_i0.astype(int)
# Townsend deprivation index - continuous
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p22189.isna(), 'p22189'] = np.nanmedian(
    proteomics_lifestyle_eur['p22189'])
# education - categorical
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur['p6138_i0'].isna(), 'p6138_i0'] = 'no_higher_education'
proteomics_lifestyle_eur.loc[
    proteomics_lifestyle_eur['p6138_i0'].str.contains('College or University degree'), 'p6138_i0'] = 'higher_education'
proteomics_lifestyle_eur.loc[
    ~proteomics_lifestyle_eur['p6138_i0'].str.contains('higher_education'), 'p6138_i0'] = 'no_higher_education'
# physical activity - continuous
proteomics_lifestyle_eur.loc[~proteomics_lifestyle_eur['p884_i0'].isin([str(x) for x in range(8)]), 'p884_i0'] = round(
    np.nanmedian(proteomics_lifestyle_eur[proteomics_lifestyle_eur['p884_i0'].isin([str(x) for x in range(8)])][
                     'p884_i0'].astype(int)))
proteomics_lifestyle_eur.loc[:, 'p884_i0'] = proteomics_lifestyle_eur.p884_i0.astype(int)
# grip strength - continuous
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p46_i0.isna(), 'p46_i0'] = np.nanmedian(
    proteomics_lifestyle_eur.p46_i0)
proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur.p47_i0.isna(), 'p47_i0'] = np.nanmedian(
    proteomics_lifestyle_eur.p47_i0)
proteomics_lifestyle_eur['grip_strength'] = (proteomics_lifestyle_eur.p46_i0 + proteomics_lifestyle_eur.p47_i0) / 2
# various food items - continuous
for field in ['p1289_i0', 'p1299_i0', 'p1309_i0', 'p1319_i0', 'p1438_i0', 'p1458_i0', 'p1528_i0']:
    proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur[field].isna(), field] = 'Do not know'
    proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur[field] == 'Less than one', field] = '0.5'
    proteomics_lifestyle_eur.loc[
        proteomics_lifestyle_eur[field].isin(['Prefer not to answer', 'Do not know']), field] = round(
        np.nanmedian(proteomics_lifestyle_eur[proteomics_lifestyle_eur[field].isin([str(x) for x in range(101)])][
                         field].astype(float)))
    proteomics_lifestyle_eur.loc[:, field] = proteomics_lifestyle_eur[field].astype(float)
# various food items - categorical
for field in ['p' + str(1329 + 10 * x) + '_i0' for x in range(7)]:
    proteomics_lifestyle_eur.loc[proteomics_lifestyle_eur[field].isna(), field] = 'Do not know'
    proteomics_lifestyle_eur.loc[
        proteomics_lifestyle_eur[field].isin(['Prefer not to answer', 'Do not know']), field] = 'Unanswered'
# export
proteomics_lifestyle_eur.drop(['p48_i0', 'p49_i0'], axis=1).to_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv', index=False)

""" Filter phenotype data to new eid set """
proteomics_phenotypes = read_csv('~/data/external/proteomics/proteomics_phenotype_file_off_dnanexus.csv',
                                 low_memory=True)
proteomics_phenotypes = proteomics_phenotypes[proteomics_phenotypes['eid'].isin(list(covars['eid']))]
proteomics_phenotypes.to_csv('~/data/internal/proteomics/proteomics/proteomics_pan_ukbb_eur_phenotypes.csv',
                             index=False)
