""" Filter phenotype data obtained off DNAnexus to get all relevant phenotype files for a given omic. """
from pandas import read_csv, concat, merge
from datetime import datetime
from pyprind import ProgBar

import sys
import argparse


# bash argument parser
parser = argparse.ArgumentParser(description='Which phenotypes to extract.')
parser.add_argument('-f', '--phenotype_dnanexus_file', type=str, required=True, help='Path to DNAnexus phenotype data to be filtered.')
parser.add_argument('-c', '--phenotype_covariate_file', type=str, required=False, default=None, help='Path to CSV covariate file containing EIDS and date attending the assessment centre (p53_i0). Usually the technical covariate file.')
parser.add_argument('-k', '--eids_to_keep', type=str, required=False, default=None, help='Path to one-column TSV with no header listing EIDs to keep. Used when no covariate file is supplied.')
parser.add_argument('-r', '--eids_to_remove', type=str, required=False, default=None, help='Path to one-column TSV with no header listing EIDs to remove. Used when no covariate file is supplied.')
parser.add_argument('-e', '--export_eids_filtered', required=False, action='store_true', help='Provide to export the phenotype file only filtered to desired EIDs. If provided, no other operations executed.')
parser.add_argument('-s', '--phenotype_set', type=str, required=True, help='Path to phenotype set, e.g. common_unisex, ard. Single column of date fields with "disease_fields" as header, in CSV format.')
parser.add_argument('-t', '--phenotype_time', type=str, required=True, choices=['before_blood', 'before_disease', 'all'], help='Phenotype time set: before_blood (those diagnosed before assessment centre), before_disease (those diagnosed before a given disease), all.')
parser.add_argument('-o', '--output', type=str, required=True, help='Destination path for exported file.')
args = parser.parse_args()  # parse arguments
phenotype_file = args.phenotype_dnanexus_file
covariate_file = args.phenotype_covariate_file
eids_to_keep_file = args.eids_to_keep
eids_to_remove_file = args.eids_to_remove
export_eids_filtered = args.export_eids_filtered
phenotype_set = args.phenotype_set
phenotype_time_set = args.phenotype_time
output_file = args.output
print(f"Options: \n\t- covariate file: {covariate_file} "
      f"\n\t- EIDs to keep: {eids_to_keep_file} "
      f"\n\t- EIDs to remove: {eids_to_remove_file} "
      f"\n\t- exporting file with only EIDS filtered: {export_eids_filtered} "
      f"\n\t- phenotype set file: {phenotype_set} "
      f"\n\t- phenotype time set: {phenotype_time_set} "
      f"\n\t- output: {output_file}",
      flush=True)

# read in data
# metabolomics_phenotypes = read_csv('~/metabolomics/metabolomics_phenotype_file_off_dnanexus.csv', low_memory=True)
phenotypes = read_csv(phenotype_file, low_memory=False)
disease_set = read_csv(phenotype_set)
# qu10_50_diseases = read_csv('~/pan_ukbb_eur_qu10_50_disease_fields.txt')
# common_unisex_date_fields = read_csv('~/common_unisex_disease_date_fields.csv')

# filter from file downloaded off DNAnexus to get Pan-UKBB EUR phenotypes
if covariate_file is not None:
    # metabolomics_covars = read_csv('~/metabolomics/metabolomics_pan_ukbb_eur_covars.csv')
    covars = read_csv(covariate_file, low_memory=False)
    phenotypes = phenotypes[phenotypes['eid'].isin(list(covars['eid']))]
else:
    if eids_to_keep_file is not None:
        eids_to_keep = read_csv(eids_to_keep_file, header=None, sep='\t')[0].to_list()
        phenotypes = phenotypes[phenotypes['eid'].isin(eids_to_keep)].copy()
    else:
        eids_to_keep = []
    if eids_to_remove_file is not None:
        eids_to_remove = read_csv(eids_to_remove_file, header=None, sep='\t')[0].to_list()
        phenotypes = phenotypes[~phenotypes['eid'].isin(eids_to_remove)].copy()
    else:
        eids_to_remove = []
    if len(eids_to_remove) * len(eids_to_keep) == 0:
        print(f"Warning: no EIDs filtering has been performed due to no covariate or EIDs files provided.", flush=True)

if export_eids_filtered == 'y':
    phenotypes.to_csv(output_file, index=False)
    print(f"{phenotypes.shape[1] - 1} phenotypes for {phenotypes.shape[0]} samples exported to {output_file}.",
          flush=True)
    sys.exit()

# filter to get common-unisex date fields
phenotypes_filtered = phenotypes[['eid'] + list(disease_set.disease_field)]

if phenotype_time_set == 'all':
    phenotypes_filtered_no_eid = phenotypes_filtered.drop('eid', axis=1)
    phenotypes_filtered_no_eid = phenotypes_filtered_no_eid.astype('O')
    phenotypes_filtered_no_eid[phenotypes_filtered_no_eid.isna()] = '0'
    phenotypes_filtered_no_eid[phenotypes_filtered_no_eid != '0'] = '1'
    phenotypes_filtered_no_eid = phenotypes_filtered_no_eid.astype('int8')
    phenotypes_filtered = concat([phenotypes_filtered[['eid']],
                                               phenotypes_filtered_no_eid], axis=1)

    phenotypes_filtered.to_csv(output_file, index=False)

    print(f"{phenotypes_filtered.shape[1] - 1} phenotypes for {phenotypes_filtered.shape[0]} samples exported to {output_file}.",
          flush=True)
    sys.exit()

elif phenotype_time_set == 'before_blood':
    if covariate_file is None:
        raise ValueError(f"A phenotype_covariate_file is required if prior phenotypes are to be obtained. Exiting.")
    # add date of attending the assessment centre
    data_df = merge(phenotypes_filtered, covars[['eid', 'p53_i0']], on='eid').reset_index(
        drop=True)
    assert (len(phenotypes_filtered) == len(data_df))

    prior_data_df = data_df[['eid']].copy()
    bar = ProgBar(iterations=phenotypes_filtered.shape[1] - 1, stream=sys.stdout, title='Processing')
    for disease in list(phenotypes_filtered.columns)[1:]:
        disease_df = data_df[data_df[disease].notnull()][['eid', 'p53_i0', disease]].reset_index(drop=True)
        for i in range(len(disease_df)):
            if any(x in disease_df.loc[i, disease] for x in ['1900-01-01', '1901-01-01', '2037-07-07',
                                                                         'Code has no event date',
                                                                         'Code has event date before participant\'s date of birth',
                                                                         'Code has event date in the future and is presumed to be a place-holder or other system default']):
                disease_df.loc[i, disease] = '0'  # phenotype is set to '0' given any of the above 'error' codes
            elif any(x in disease_df.loc[i, disease] for x in ['1902-02-02', '1903-03-03',
                                                                           'Code has event date matching participant\'s date of birth',
                                                                           'Code has event date after participant\'s date of birth and falls in the same calendar year as date of birth']):
                disease_df.loc[i, disease] = '1'  # phenotype is set to '1' given any of the above 'error' codes
            else:
                date_of_assessment = datetime.strptime(disease_df.loc[i, 'p53_i0'], '%Y-%m-%d')
                date_of_onset = datetime.strptime(disease_df.loc[i, disease], '%Y-%m-%d')

                if date_of_onset >= date_of_assessment:  # onset is on or after assessment
                    disease_df.loc[i, disease] = '0'  # phenotype is removed (set to '0')

        prior_data_df = merge(prior_data_df, disease_df[['eid', disease]], on='eid', how='left')
        bar.update()

    phenotypes_only = prior_data_df.iloc[:, 1:].copy()
    phenotypes_only = phenotypes_only.astype('O')
    phenotypes_only[phenotypes_only.isna()] = '0'
    phenotypes_only[phenotypes_only != '0'] = '1'
    phenotypes_only = phenotypes_only.astype('int8')

    final_df = concat([prior_data_df['eid'], phenotypes_only], axis=1)
    final_df.to_csv(output_file, index=False)

    print(f"{final_df.shape[1] - 1} phenotypes for {final_df.shape[0]} samples exported to {output_file}.", flush=True)
    sys.exit()

elif phenotype_time_set == 'before_disease':
    pass
else:
    print(f"Phenotype time set {phenotype_time_set} is not recognised. Should be one of ['prior', 'all'].", flush=True)
