""" General script to be submitted via bash to generate principal components representing ~50% variance explained by
prior/accompanying/comorbid (however you like to describe it) disease diagnoses.
Allows for disease diagnoses to be resolved by time relative to events, e.g. blood draw, diagnosis of an ARD. """
from pandas import read_csv, concat, merge, to_datetime, DataFrame
from sklearn.decomposition import PCA
from pyprind import ProgBar

import sys
import argparse
import seaborn as sns
import pickle
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# bash argument parser
parser = argparse.ArgumentParser(description='Which phenotypes to extract.')
parser.add_argument('--phenotype_dnanexus_file', type=str, required=True,
                    help='Path to DNAnexus phenotype data to be filtered.')
parser.add_argument('--phenotype_covariate_file', type=str, required=False, default=None,
                    help='Path to CSV covariate file containing EIDS and date attending the assessment centre (p53_i0). Usually the technical covariate file.')
parser.add_argument('--eids_to_keep', type=str, required=False, default=None,
                    help='Path to one-column TSV with no header listing EIDs to keep. Used when no covariate file is supplied.')
parser.add_argument('--eids_to_remove', type=str, required=False, default=None,
                    help='Path to one-column TSV with no header listing EIDs to remove. Used when no covariate file is supplied.')
# parser.add_argument('-s', '--phenotype_set', type=str, required=True, help='Path to phenotype set, e.g. common_unisex, ard. Single column of date fields with "disease_fields" as header, in CSV format.')
parser.add_argument('--background_phenotype_set', type=str, required=True,
                    help='Path to set of phenotypes from which PCs are to be derived. One-column CSV.')
parser.add_argument('--diagnosis_time_relation_to_ards', required=False, default=None,
                    choices=['before', 'after'],
                    help='Whether diagnoses used to generated PCs should be "before" or "after" ARD diagnoses. Default = all diagnoses used.')
parser.add_argument('--diagnosis_time_relation_to_blood_draw', required=False, default=None,
                    choices=['before', 'after'],
                    help='Whether diagnoses used to generate PCs should be "before" or "after" blood draw. Default = all diagnoses used.')
parser.add_argument('--pca_export_dir', type=str, required=True,
                    help='Destination path for exported principal components.')
parser.add_argument('--plot_clustermap', required=False, action='store_true',
                    help='Include flag to plot a clustermap clustering features\' weights on components should be plotted. Requires --clustermap_path.')
parser.add_argument('--clustermap_path', type=str, required=False, default=None,
                    help='Destination path for exported clustermap plot, without file extension.')
parser.add_argument('--cumulative_variance_path', type=str, required=False, default=None,
                    help='Destination path for selected PCs\' cumulative variance explained plot, without file extension.')
args = parser.parse_args()  # parse arguments
phenotype_file = args.phenotype_dnanexus_file
covariate_file = args.phenotype_covariate_file
eids_to_keep_file = args.eids_to_keep
eids_to_remove_file = args.eids_to_remove
background_phenotype_set = args.background_phenotype_set
diagnosis_time_relation_to_ards = args.diagnosis_time_relation_to_ards
diagnosis_time_relation_to_blood_draw = args.diagnosis_time_relation_to_blood_draw
pca_export_dir = args.pca_export_dir
plot_clustermap = args.plot_clustermap
clustermap_path = args.clustermap_path
cumulative_variance_path = args.cumulative_variance_path
print(f"Options: "
      f"\n\t- phenotype file: {phenotype_file} "
      f"\n\t- covariate file: {covariate_file} "
      f"\n\t- EIDs to keep: {eids_to_keep_file} "
      f"\n\t- EIDs to remove: {eids_to_remove_file} "
      f"\n\t- background phenotype set: {background_phenotype_set} "
      f"\n\t- diagnosis time relation to ARDs: {diagnosis_time_relation_to_ards} "
      f"\n\t- diagnosis time relation to blood draw: {diagnosis_time_relation_to_blood_draw} "
      f"\n\t- PCs directory: {pca_export_dir}"
      f"\n\t- whether to plot a clustermap of PCs: {plot_clustermap} "
      f"\n\t- clustermap export directory: {clustermap_path} "
      f"\n\t- plot of selected PCs' cumulative variance explained: {cumulative_variance_path}",
      flush=True)

# some data input validation
if diagnosis_time_relation_to_blood_draw is not None:
    assert covariate_file is not None, 'A covariate file containing dates of assessment centre visit is required to get diagnoses before or after blood draw.'
if plot_clustermap:
    assert clustermap_path is not None, 'Destination path for exported clustermap plot is required if a clustermap is to be plotted.'

# read in phenotype file
phenotypes = read_csv(phenotype_file, low_memory=False)

# read in background disease set - the diseases whose diagnoses will be collected (usually common-unisex diseases)
disease_set = read_csv(background_phenotype_set)

# read in ARD diseases
qu10_50_diseases = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')

# filter EIDs
if covariate_file is not None:
    covars = read_csv(covariate_file, low_memory=False)
    phenotypes = phenotypes[phenotypes['eid'].isin(list(covars['eid']))].reset_index(drop=True)
else:
    if eids_to_keep_file is not None:
        eids_to_keep = read_csv(eids_to_keep_file, header=None, sep='\t')[0].to_list()
        phenotypes = phenotypes[phenotypes['eid'].isin(eids_to_keep)].reset_index(drop=True).copy()
    else:
        eids_to_keep = []
    if eids_to_remove_file is not None:
        eids_to_remove = read_csv(eids_to_remove_file, header=None, sep='\t')[0].to_list()
        phenotypes = phenotypes[~phenotypes['eid'].isin(eids_to_remove)].reset_index(drop=True).copy()
    else:
        eids_to_remove = []
    if len(eids_to_remove) * len(eids_to_keep) == 0:
        print(f"Warning: no EIDs filtering has been performed due to no covariate or EIDs files provided.", flush=True)
# filter to get common-unisex date fields
phenotypes_filtered = phenotypes[['eid'] + list(disease_set.disease_field)]

# replace error codes to facilitate diagnosis time comparisons
# these are taken to be non-diagnoses, so are set to a date to be turned to 0 later on
codes_for_non_diagnosis = ['1900-01-01', '1901-01-01', '2037-07-07', 'Code has no event date',
                           'Code has event date before participant\'s date of birth',
                           'Code has event date in the future and is presumed to be a place-holder or other system default']
codes_for_non_diagnosis_pattern = '|'.join(codes_for_non_diagnosis)
# these are taken to be early-in-life diagnoses, so are set to a date to be turned to 1 later on
codes_for_early_diagnosis = ['1902-02-02', '1903-03-03',
                             'Code has event date matching participant\'s date of birth',
                             'Code has event date after participant\'s date of birth and falls in the same calendar year as date '
                             'of birth']
codes_for_early_diagnosis_pattern = '|'.join(codes_for_early_diagnosis)

# perform disease-agnostic phenotype filtering
phenotypes_to_pc = phenotypes_filtered.copy()
# if wanting diseases before blood draw
if diagnosis_time_relation_to_blood_draw is not None:
    data_df = merge(phenotypes_filtered, covars[['eid', 'p53_i0']], on='eid')
    if diagnosis_time_relation_to_blood_draw == 'before':
        # print('blood-draw-before', flush=True)
        phenotypes_dates_no_eid = (data_df.drop(columns='eid')
                                   .apply(
            # want to remove these dates, so set them as 'after'
            lambda col: col.mask(col.str.contains(codes_for_non_diagnosis_pattern, na=False), '2037-07-07'))
                                   .apply(
            # want to keep these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_early_diagnosis_pattern, na=False), '1902-02-02'))
                                   # empty cells are taken to be non-diagnoses, so are set to a date in the future to be turned to 0 later on
                                   .fillna('2037-07-07')
                                   # convert phenotypes and assessment centre dates to datetime format
                                   .apply(to_datetime))
        # get phenotypes diagnosed before assessment centre date
        phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns='p53_i0').lt(phenotypes_dates_no_eid['p53_i0'],
                                                                             axis=0)
    elif diagnosis_time_relation_to_blood_draw == 'after':
        # print('blood-draw-after', flush=True)
        phenotypes_dates_no_eid = (data_df.drop(columns='eid')
                                   .apply(
            # want to remove these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_non_diagnosis_pattern, na=False), '1902-02-02'))
                                   .apply(
            # want to remove these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_early_diagnosis_pattern, na=False), '1902-02-02'))
                                   # empty cells are taken to be non-diagnoses, so are set to a date in the past to be turned to 0 later on
                                   .fillna('1902-02-02')
                                   # convert phenotypes and assessment centre dates to datetime format
                                   .apply(to_datetime))
        phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns='p53_i0').gt(phenotypes_dates_no_eid['p53_i0'],
                                                                             axis=0)
    # add back EID column
    # phenotypes_to_pc = concat([phenotypes_filtered[['eid']], phenotypes_to_pc])
    print(phenotypes_to_pc.sum().sum(), flush=True)

# set up phenotype data for disease-specific phenotype filtering
# if wanting diseases before/after ARD diagnoses
if diagnosis_time_relation_to_ards is not None:
    if diagnosis_time_relation_to_blood_draw is None:
        phenotypes_to_pc.drop(columns='eid', inplace=True)
    else:
        # if time comparison to blood draw is done, then phenotypes_to_pc is a boolean mask
        # obtain diagnosis dates that fulfill the comparison to the blood draw date
        phenotypes_to_pc = phenotypes_filtered.drop(columns='eid')[phenotypes_to_pc]

    if diagnosis_time_relation_to_ards == 'before':
        # print('ard-before', flush=True)
        phenotypes_dates_no_eid = (phenotypes_to_pc
                                   .apply(
            # want to remove these dates, so set them as 'after'
            lambda col: col.mask(col.str.contains(codes_for_non_diagnosis_pattern, na=False), '2037-07-07'))
                                   .apply(
            # want to keep these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_early_diagnosis_pattern, na=False), '1902-02-02'))
                                   # empty cells are taken to be non-diagnoses, so are set to a date in the future to be turned to 0 later on
                                   .fillna('2037-07-07')
                                   # convert phenotypes to datetime format
                                   .apply(to_datetime))
    elif diagnosis_time_relation_to_ards == 'after':
        # print('ard-after', flush=True)
        phenotypes_dates_no_eid = (phenotypes_to_pc
                                   .apply(
            # want to remove these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_non_diagnosis_pattern, na=False), '1902-02-02'))
                                   .apply(
            # want to remove these dates, so set them as 'before'
            lambda col: col.mask(col.str.contains(codes_for_early_diagnosis_pattern, na=False), '1902-02-02'))
                                   # empty cells are taken to be non-diagnoses, so are set to a date in the past to be turned to 0 later on
                                   .fillna('1902-02-02')
                                   # convert phenotypes to datetime format
                                   .apply(to_datetime))

# if no diagnosis_time_relation has been set, then phenotypes_to_pc is just phenotypes_filtered
# convert into a boolean mask, where valid diagnoses are True and everything else is False
if diagnosis_time_relation_to_blood_draw is None and diagnosis_time_relation_to_ards is None:
    phenotypes_to_pc = (phenotypes_filtered.drop(columns='eid')
                        .apply(
        # want to remove these dates, so set them to nan
        lambda col: col.mask(col.str.contains(codes_for_non_diagnosis_pattern, na=False), np.nan))
                        .apply(
        # want to remove these dates, so set them to nan
        lambda col: col.mask(col.str.contains(codes_for_early_diagnosis_pattern, na=False), np.nan)))
    phenotypes_to_pc = phenotypes_to_pc.notna()

# df showing cumulative count of %variance explained by first 20 PCs for each ARD
all_var_exp_df = DataFrame()
components_across_disease = []
common_unisex_diseases_list = read_csv('~/data/internal/phenotype_coding/common_unisex_disease_date_fields.csv')['disease_field'].to_list()
bar = ProgBar(len(qu10_50_diseases.disease_field), stream=sys.stdout, title='Getting previous diagnoses and PCA-ing')
for disease in qu10_50_diseases.disease_field:

    # if wanting diseases before/after ARD diagnoses
    if diagnosis_time_relation_to_ards is not None:

        if diagnosis_time_relation_to_ards == 'before':
            phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns=disease).lt(phenotypes_dates_no_eid[disease],
                                                                                axis=0)
            # account for non-diagnoses of ARDs, which will have been converted to '2037-07-07' in the ARD column and ellicited lots of True in the dataframe
            # all True in those rows should be False because the ARD itself has not been diagnosed, so cannot have disease diagnosed before it
            phenotypes_to_pc.loc[phenotypes_dates_no_eid[disease] == '2037-07-07', :] = False
        elif diagnosis_time_relation_to_ards == 'after':
            phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns=disease).gt(phenotypes_dates_no_eid[disease],
                                                                                axis=0)
            # non-diagnoses of ARDs, which will have been converted to '2037-07-07' in the ARD column will ellicit lots of 0s in the dataframe
            phenotypes_to_pc.loc[phenotypes_dates_no_eid[disease] == '1902-02-02', :] = False

        print(phenotypes_to_pc.sum().sum(), flush=True)

    # phenotypes_to_pc should be just the phenotypes in mask form, with no 'eid' or 'p53_i0' column
    if 'eid' in phenotypes_to_pc.columns:
        phenotypes_to_pc.drop(columns='eid', inplace=True)

    # phenotypes_to_pc should not have the disease being investigated
    if disease in phenotypes_to_pc.columns:
        phenotypes_to_pc.drop(columns=disease, inplace=True)

    # one-hot the diagnoses
    phenotypes_to_pcs_one_hot = phenotypes_to_pc.astype(int)

    # fit PCA and dimension-reduce previous diagnosis binary dataset
    pca = PCA(random_state=42)
    phenotypes_reduced = pca.fit_transform(phenotypes_to_pcs_one_hot)

    var_exp = pca.explained_variance_ratio_
    all_var_exp_df = concat(
        [all_var_exp_df, concat([DataFrame({'disease': [disease]}), DataFrame([var_exp])], axis=1)])

    bar.update()

all_var_exp_df.columns = ['disease'] + ['PC' + str(x) for x in list(range(1, 306))]
# turn all_var_exp_df into a cumulative dataset
all_var_exp_cum_df = all_var_exp_df.copy()
for pc in range(2, 306):
    all_var_exp_cum_df['PC' + str(pc)] = all_var_exp_cum_df['PC' + str(pc)] + all_var_exp_cum_df['PC' + str(pc - 1)]

# choose how many PCs to include based on the maximum number it takes to accrue 50% variance explained
selected_var_exp_cum_df = all_var_exp_cum_df.drop('disease', axis=1).loc[:,
                          (all_var_exp_cum_df.drop('disease', axis=1) < 0.5).any()]
num_pcs_selected = selected_var_exp_cum_df.shape[1]
avg_var_explained_by_selected_pcs = all_var_exp_cum_df[f"PC{num_pcs_selected}"].mean()
print(all_var_exp_cum_df[f"PC{num_pcs_selected}"].describe())

# plot cumulative variance explained by the selected number of PCs
disease_info = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
plt.figure(figsize=(8, 6))
sns.barplot(data=all_var_exp_cum_df[['disease', f"PC{num_pcs_selected}"]].rename(columns={f"PC{num_pcs_selected}": 'Fraction of variance explained'})
            .merge(disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'})),
            x='code_chapter', y='Fraction of variance explained')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.gca().tick_params(axis='x', rotation=90, labelsize=8)
plt.tight_layout()
plt.savefig(f"{cumulative_variance_path}.svg")
plt.savefig(f"{cumulative_variance_path}.png")

# components_across_disease = []
# common_unisex_diseases_list = read_csv('~/common_unisex_disease_date_fields.csv')['disease_field'].to_list()
bar = ProgBar(len(qu10_50_diseases.disease_field), stream=sys.stdout, title='Exporting relevant PCA results')
for i, disease in enumerate(qu10_50_diseases.disease_field):

    # if wanting diseases before/after ARD diagnoses
    if diagnosis_time_relation_to_ards is not None:
        if diagnosis_time_relation_to_ards == 'before':
            phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns=disease).lt(phenotypes_dates_no_eid[disease],
                                                                                axis=0)
            # account for non-diagnoses of ARDs, which will have been converted to '2037-07-07' in the ARD column and ellicited lots of True in the dataframe
            # all True in those rows should be False because the ARD itself has not been diagnosed, so cannot have disease diagnosed before it
            phenotypes_to_pc.loc[phenotypes_dates_no_eid[disease] == '2037-07-07', :] = False
        elif diagnosis_time_relation_to_ards == 'after':
            phenotypes_to_pc = phenotypes_dates_no_eid.drop(columns=disease).gt(phenotypes_dates_no_eid[disease],
                                                                                axis=0)
            # non-diagnoses of ARDs, which will have been converted to '2037-07-07' in the ARD column will ellicit lots of 0s in the dataframe
            phenotypes_to_pc.loc[phenotypes_dates_no_eid[disease] == '1902-02-02', :] = False

    # phenotypes_to_pc should be just the phenotypes in mask form, with no 'eid' or 'p53_i0' column
    if 'eid' in phenotypes_to_pc.columns:
        phenotypes_to_pc.drop(columns='eid', inplace=True)

    # phenotypes_to_pc should not have the disease being investigated
    if disease in phenotypes_to_pc.columns:
        phenotypes_to_pc.drop(columns=disease, inplace=True)

    # one-hot the diagnoses
    phenotypes_to_pcs_one_hot = phenotypes_to_pc.astype(int)

    # fit PCA and dimension-reduce previous diagnosis binary dataset
    pca = PCA(random_state=42, n_components=num_pcs_selected)
    phenotypes_reduced = pca.fit_transform(phenotypes_to_pcs_one_hot)

    # export if a destination is specified
    if pca_export_dir is not None:
        os.makedirs(pca_export_dir, exist_ok=True)
        (concat([phenotypes_filtered[['eid']], DataFrame(phenotypes_reduced)], axis=1)
         .set_axis(['eid'] + ['PC' + str(x) for x in range(1, num_pcs_selected + 1)], axis=1)
         .to_csv(pca_export_dir + '/' + disease + '.csv', index=False))

    components_across_disease.append(
        # inserts a column of 0s for the ARD excluded in this PCA
        np.insert(pca.components_, obj=common_unisex_diseases_list.index(disease), values=0, axis=1))

    bar.update()
