"""
Regress out certain covariate and lifestyle factors from the proteomics data. Done in prep for Cox regressions.
Assess the effect each regression has on variance explained by the lifestyle factors.
"""
from pandas import read_csv, merge, concat, DataFrame
from datetime import datetime
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, robust_scale
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

import sys

input_data = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_data.csv')

participant_covariates = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv').drop('p23099_i0', axis=1)
lifestyles = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv')
participant_covariates = merge(participant_covariates, lifestyles, on='eid')

protein_covariates = read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv')

# participant covars to regress out: well, plate row, plate column, all fish, all meat bar poultry, water
participant_covariates_of_interest = ['plate_row', 'plate_column', 'p30902_i0', 'p1329_i0', 'p1339_i0', 'p1349_i0',
                                      'p1369_i0', 'p1379_i0', 'p1389_i0', 'p1528_i0']
# plate ID must be included to facilitate merging with protein-based covariates later
# date required to get storage_time
participant_covariates = participant_covariates[['eid', 'p53_i0', 'p30901_i0'] + participant_covariates_of_interest].copy()

# get protein fields
protein_fields = list(input_data.columns)[1:]
protein = protein_fields[
    int(sys.argv[1]) - 1]  # use passed-in SGE_TASK_ID to isolate out of which protein to regress effects


def regress_out(covariate_df, data_df, covariate, protein, covariate_type='categorical'):
    """
    Regress out the effect of a covariate from a protein's data.
    Also calculate the effect of that covariate on that protein's variance before and after regressing it out.
    If the covariate is categorical, it is transformed by one-hot encoding before regression.
    :param covariate_df:
    :param data_df:
    :param covariate:
    :param protein:
    :param covariate_type:
    :return:
    """
    covariate_data_df = merge(covariate_df[['eid', covariate]],
                              data_df[['eid', protein]], on='eid')

    covariate_data_df = covariate_data_df[covariate_data_df[protein].notna()].reset_index(drop=True)

    if covariate_type == 'categorical':
        # define the column transformer for one-hot encoding
        categorical_features = [covariate]
        preprocessor = ColumnTransformer(
            transformers=[('cat', OneHotEncoder(), categorical_features)],
            remainder='passthrough'
        )

        # create the pipeline, incorporating the pre-processing step and a linear regressor
        pipeline_to_regress_out = Pipeline(steps=[('preprocessor', preprocessor), ('model', LinearRegression())])
        pipeline_to_check_after_regression = Pipeline(
            steps=[('preprocessor', preprocessor), ('model', LinearRegression())])

    elif covariate_type == 'continuous':
        pipeline_to_regress_out = LinearRegression()
        pipeline_to_check_after_regression = LinearRegression()

    else:
        raise ValueError('covariate_type must be one of \'categorical\' and \'continuous\'')

    # isolate covariate and protein
    x = covariate_data_df[[covariate]]
    y = covariate_data_df[[protein]]

    # fit the pipeline
    pipeline_to_regress_out.fit(x, y)

    # score the fit
    var_explained_before = pipeline_to_regress_out.score(x, y)

    # obtain predictions
    y_pred = pipeline_to_regress_out.predict(x)

    # obtain residuals
    residuals = y - y_pred

    # fit pipeline again and check variance explained
    var_explained_after = pipeline_to_check_after_regression.fit(x, residuals).score(x, residuals)

    return var_explained_before, concat([covariate_data_df[['eid']], residuals], axis=1), var_explained_after


def var_explained_by_one_covar(covars_df, response_var_df, covar, response_var, covar_type='categorical'):
    """
    Determine variance of the response variable explained by the technical covariate through linear regression.
    If the technical covariate is categorical, it is transformed by one-hot encoding before regression.
    :param covars_df:
    :param response_var_df:
    :param covar:
    :param response_var:
    :param covar_type:
    :return:
    """
    covar_response_var = merge(covars_df[['eid', covar]],
                               response_var_df, on='eid')

    covar_response_var = covar_response_var[covar_response_var[response_var].notna()].reset_index(drop=True)

    if covar_type == 'categorical':
        # define the column transformer for one-hot encoding
        categorical_features = [covar]
        preprocessor = ColumnTransformer(
            transformers=[('cat', OneHotEncoder(), categorical_features)],
            remainder='passthrough'
        )

        # create the pipeline, incorporating the pre-processing step and a linear regressor
        pipeline = Pipeline(steps=[('preprocessor', preprocessor), ('model', LinearRegression())])

    elif covar_type == 'continuous':
        pipeline = LinearRegression()

    else:
        raise ValueError('covar_type must be one of \'categorical\' and \'continuous\'')

    # isolate explanatory and response variables
    x = covar_response_var[[covar]]
    y = covar_response_var[[response_var]]

    # fit the pipeline and score the fit
    var_explained = pipeline.fit(x, y).score(x, y)

    return var_explained


# covariates should be of form [{'field': <field>, 'name': <name>, 'type': 'continuous'/'categorical'}, {}, ...]
def var_explained_by_each_covar(covariate_df, data_df, covariates, protein, stage):
    vars_explained_by_covars = {'stage': [stage] * len(covariates),
                                'covar_name': [],
                                'var_explained': [],
                                'protein': [protein] * len(covariates)}
    for covar_dict in covariates:
        vars_explained_by_covars['covar_name'].append(covar_dict['name'])
        vars_explained_by_covars['var_explained'].append(var_explained_by_one_covar(covars_df=covariate_df,
                                                                                    response_var_df=data_df,
                                                                                    covar=covar_dict['field'],
                                                                                    response_var=protein,
                                                                                    covar_type=covar_dict['type']))
    return DataFrame.from_dict(vars_explained_by_covars)


covariate_info = [{'field': 'p30902_i0', 'name': 'Well', 'type': 'categorical'},
                  {'field': 'Processing_StartDate', 'name': 'Date processed', 'type': 'categorical'},
                  {'field': 'storage_time', 'name': 'Storage time', 'type': 'continuous'},
                  {'field': 'plate_row', 'name': 'Plate row', 'type': 'categorical'},
                  {'field': 'plate_column', 'name': 'Plate column', 'type': 'categorical'},
                  {'field': 'p1329_i0', 'name': 'Oily fish intake', 'type': 'categorical'},
                  {'field': 'p1339_i0', 'name': 'Non-oily fish intake', 'type': 'categorical'},
                  {'field': 'p1349_i0', 'name': 'Processed meat intake', 'type': 'categorical'},
                  {'field': 'p1369_i0', 'name': 'Beef intake', 'type': 'categorical'},
                  {'field': 'p1379_i0', 'name': 'Lamb intake', 'type': 'categorical'},
                  {'field': 'p1389_i0', 'name': 'Pork intake', 'type': 'categorical'},
                  {'field': 'p1528_i0', 'name': 'Water intake', 'type': 'continuous'}]
vars_explained_by_all_covars_at_each_step = DataFrame(columns=['stage', 'covar_name', 'var_explained', 'protein'])

if protein in ['ervv_1', 'hla_a', 'hla_dra', 'hla_e']:  # these proteins have '-' instead of '_' in the protein covars
    single_protein_covars = protein_covariates[
        protein_covariates['Assay'].str.lower() == protein.replace('_', '-')].copy()
else:
    single_protein_covars = protein_covariates[
        protein_covariates['Assay'].str.lower() == protein].copy()

eid_plate_protein_covars = merge(participant_covariates[['eid', 'p30901_i0']],
                                 single_protein_covars.drop(['Assay', 'UniProt'], axis=1),
                                 left_on='p30901_i0', right_on='PlateID')

all_covars = merge(participant_covariates,
                   eid_plate_protein_covars.drop(['p30901_i0', 'PlateID', 'Panel', 'Panel_Lot_Nr'], axis=1),
                   on='eid')

# storage_time not something to regress out
# calculate time between sample taken and sample processed
times_between_sample_taken_and_measurement = []
for i in range(len(all_covars)):
    times_between_sample_taken_and_measurement.append((datetime.strptime(
        all_covars.loc[i, 'Processing_StartDate'], '%Y-%m-%d') - datetime.strptime(
        all_covars.loc[i, 'p53_i0'], '%Y-%m-%d')).days)
all_covars['storage_time'] = times_between_sample_taken_and_measurement

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=input_data,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='before (0)')])

# well - step 1
step1_var_explained_before, step1_residuals, step1_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=input_data,
                                                                                     covariate='p30902_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step1_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='1')])

# date at centre - step 2
step2_var_explained_before, step2_residuals, step2_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step1_residuals,
                                                                                     covariate='plate_row',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step2_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='2')])

# plate - step 3
step3_var_explained_before, step3_residuals, step3_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step2_residuals,
                                                                                     covariate='plate_column',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step3_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='3')])

# processing date - step 4
step4_var_explained_before, step4_residuals, step4_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step3_residuals,
                                                                                     covariate='Processing_StartDate',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step4_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='4')])

# oily fish - step 5
step5_var_explained_before, step5_residuals, step5_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step4_residuals,
                                                                                     covariate='p1329_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step5_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='5')])

# non-oily fish - step 6
step6_var_explained_before, step6_residuals, step6_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step5_residuals,
                                                                                     covariate='p1339_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step6_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='6')])

# processed meat - step 7
step7_var_explained_before, step7_residuals, step7_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step6_residuals,
                                                                                     covariate='p1349_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step7_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='7')])

# beef - step 8
step8_var_explained_before, step8_residuals, step8_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step7_residuals,
                                                                                     covariate='p1369_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step8_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='8')])

# lamb - step 9
step9_var_explained_before, step9_residuals, step9_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                     data_df=step8_residuals,
                                                                                     covariate='p1379_i0',
                                                                                     protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step9_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='9')])

# lamb - step 10
step10_var_explained_before, step10_residuals, step10_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                        data_df=step9_residuals,
                                                                                        covariate='p1389_i0',
                                                                                        protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step10_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='10')])

# water - step 11
step11_var_explained_before, step11_residuals, step11_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                        data_df=step10_residuals,
                                                                                        covariate='p1528_i0',
                                                                                        protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=step11_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='11')])

# storage time - step 12
step12_var_explained_before, final_residuals, step12_var_explained_after = regress_out(covariate_df=all_covars,
                                                                                       data_df=step11_residuals,
                                                                                       covariate='storage_time',
                                                                                       protein=protein)

vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                    var_explained_by_each_covar(covariate_df=all_covars,
                                                                                data_df=final_residuals,
                                                                                covariates=covariate_info,
                                                                                protein=protein,
                                                                                stage='after (12)')])

final_residuals.to_csv('regressing_out_for_cox/post_regression_protein_data/' + protein + '.csv', index=False)
vars_explained_by_all_covars_at_each_step.to_csv(
    'regressing_out_for_cox/effects/' + protein + '_effects_during_regressions.csv',
    index=False)
