""" Perform SMOTE-RUS and robust scaling on protein data before associating to disease with logistic regression. """
from pandas import read_csv, merge, concat, to_datetime
from imblearn.over_sampling import SMOTENC
from imblearn.under_sampling import RandomUnderSampler
from sklearn.preprocessing import robust_scale

import numpy as np
import sys

np.random.seed(42)

""" Read in command line argument (which should be a protein field only)"""
command_line_args = sys.argv
protein = command_line_args[1]

""" Read in other data """
data = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_post_regressions_data.csv')
participant_covariates = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv')
protein_covariates = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv')
phenotypes = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv')
diseases = list(read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])

# create mean date sample processed
protein_covariates['processing_start_datetime'] = to_datetime(protein_covariates['Processing_StartDate'])
protein_covars_mean_dates = protein_covariates.groupby('PlateID')['processing_start_datetime'].mean().reset_index().rename(columns={'PlateID': 'p30901_i0', 'processing_start_datetime': 'mean_processing_start_date'})

# merge into participant covariates and calculate storage time
covariates = merge(participant_covariates, protein_covars_mean_dates, on='p30901_i0')
covariates['storage_time'] = (covariates['mean_processing_start_date'] - to_datetime(covariates['p53_i0'])).dt.days

non_pca_covars_of_interest = ['p21003_i0', 'p31', 'p54_i0', 'storage_time']  # age, sex, centre, participant storage time
categorical_covariates = ['p31', 'p54_i0']

for disease in diseases:

    pca = read_csv(f"~/data/internal/genomics/pca/${disease}/${disease}_pca.eigenvec", delim_whitespace=True)
    pca = pca.drop(['#FID'] + ['PC' + str(x) for x in list(range(11, 21))], axis=1).rename(columns={'IID': 'eid'})

    single_phenotypes = phenotypes[['eid', disease]].copy()
    single_data = data[['eid', protein]].copy()

    covariates_pca = merge(covariates[['eid'] + non_pca_covars_of_interest], pca, on='eid', how='left')
    data_covariates_pca = merge(single_data, covariates_pca, on='eid')
    data_covariates_pca_phenotype = merge(data_covariates_pca, single_phenotypes, on='eid')

    data_covariates_pca_phenotype = data_covariates_pca_phenotype[~data_covariates_pca_phenotype[protein].isna()].copy()

    if 1 < data_covariates_pca_phenotype[disease].sum() < 6:
        num_neighbours = data_covariates_pca_phenotype[disease].sum()
    elif data_covariates_pca_phenotype[disease].sum() >= 6:
        num_neighbours = 5
    else:
        continue

    print(disease + ': num_neigbours = ' + str(num_neighbours))

    oversampling_sampling_frac = len(data_covariates_pca_phenotype) / ((len(data_covariates_pca_phenotype) - data_covariates_pca_phenotype[disease].sum()) * 2)
    oversample = SMOTENC(categorical_features=categorical_covariates, sampling_strategy=oversampling_sampling_frac,
                         random_state=42,
                         k_neighbors=num_neighbours - 1)
    undersample = RandomUnderSampler(sampling_strategy=1.0, random_state=42)

    x, y = oversample.fit_resample(data_covariates_pca_phenotype.drop(['eid', disease], axis=1), data_covariates_pca_phenotype[[disease]])
    x, y = undersample.fit_resample(x, y)

    balanced = concat([x, y], axis=1)
    balanced.loc[:, protein] = robust_scale(balanced.loc[:, protein])
    balanced.to_csv('~/ch3_proteomics/3.2 suitable_regression_algorithm/balanced_data/' + disease + '-' + protein + '_balanced_data.csv', index=False)
