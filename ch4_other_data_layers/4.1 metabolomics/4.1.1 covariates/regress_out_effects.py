"""
Regress out certain technical covariates from the log1p-transformed metabolomics data.
Assess the effect each regression has on variance explained by the covariates.
"""
from pandas import read_csv, merge, concat, DataFrame
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, robust_scale
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

import matplotlib.pyplot as plt
import seaborn as sns

input_data = read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_log1p_data.csv')
covariates = read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_covars.csv')

# covars to regress out: date at centre; prepped-for time; spectrometer; well; plate; sample prepared date, sample measured date
covariates_of_interest = ['p53_i0', 'prepped_for_time', 'p23650_i0', 'p23660_i0', 'p23649_i0', 'sample_prepared_date', 'sample_measured_date']
covariates = covariates[['eid'] + covariates_of_interest].copy()

# get metabolite fields
metabolite_fields = list(input_data.columns)[1:]


def regress_out(covariate_df, data_df, covariate, metabolite, covariate_type='categorical'):
    """
    Regress out the effect of a covariate from a metabolite's data.
    Also calculate the effect of that covariate on that metabolite's variance before and after regressing it out.
    If the covariate is categorical, it is transformed by one-hot encoding before regression.
    :param covariate_df:
    :param data_df:
    :param covariate:
    :param metabolite:
    :param covariate_type:
    :return:
    """
    covariate_data_df = merge(covariate_df[['eid', covariate]],
                              data_df[['eid', metabolite]], on='eid')

    covariate_data_df = covariate_data_df[covariate_data_df[metabolite].notna()].reset_index(drop=True)

    if covariate_type == 'categorical':
        # define the column transformer for one-hot encoding
        categorical_features = [covariate]
        preprocessor = ColumnTransformer(
            transformers=[('cat', OneHotEncoder(), categorical_features)],
            remainder='passthrough'
        )

        # create the pipeline, incorporating the pre-processing step and a linear regressor
        pipeline_to_regress_out = Pipeline(steps=[('preprocessor', preprocessor), ('model', LinearRegression())])
        pipeline_to_check_after_regression = Pipeline(steps=[('preprocessor', preprocessor), ('model', LinearRegression())])

    elif covariate_type == 'continuous':
        pipeline_to_regress_out = LinearRegression()
        pipeline_to_check_after_regression = LinearRegression()

    else:
        raise ValueError('covariate_type must be one of \'categorical\' and \'continuous\'')

    # isolate covariate and metabolite
    x = covariate_data_df[[covariate]]
    y = covariate_data_df[[metabolite]]

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
def var_explained_by_each_covar(covariate_df, data_df, covariates, metabolite, stage):
    vars_explained_by_covars = {'stage': [stage] * len(covariates),
                                'covar_name': [],
                                'var_explained': [],
                                'metabolite': [metabolite] * len(covariates)}
    for covar_dict in covariates:
        vars_explained_by_covars['covar_name'].append(covar_dict['name'])
        vars_explained_by_covars['var_explained'].append(var_explained_by_one_covar(covars_df=covariate_df,
                                                                                    response_var_df=data_df,
                                                                                    covar=covar_dict['field'],
                                                                                    response_var=metabolite,
                                                                                    covar_type=covar_dict['type']))
    return DataFrame.from_dict(vars_explained_by_covars)


covariate_info = [{'field': 'p23649_i0', 'name': 'plate', 'type': 'categorical'},
                  {'field': 'p23650_i0', 'name': 'spectrometer', 'type': 'categorical'},
                  {'field': 'prepped_for_time', 'name': 'prepped-for time', 'type': 'continuous'},
                  {'field': 'p53_i0', 'name': 'date at centre', 'type': 'categorical'},
                  {'field': 'sample_prepared_date', 'name': 'date sample prepared', 'type': 'categorical'},
                  {'field': 'sample_measured_date', 'name': 'date sample measured', 'type': 'categorical'},
                  {'field': 'p23660_i0', 'name': 'well', 'type': 'categorical'}]
all_final_residuals = DataFrame({'eid': list(covariates.eid)})
vars_explained_by_all_covars_at_each_step = DataFrame(columns=['stage', 'covar_name', 'var_explained', 'metabolite'])
# bar = ProgBar(iterations=len(metabolite_fields), stream=sys.stdout, title='Regressing out technical covariates')
for metabolite in metabolite_fields:
    print(metabolite, flush=True)
    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=input_data,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='before (0)')])

    # plate - step 1
    step1_var_explained_before, step1_residuals, step1_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=input_data,
                                                                                         covariate='p23649_i0',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step1_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='1')])

    # spectrometer - step 2
    step2_var_explained_before, step2_residuals, step2_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step1_residuals,
                                                                                         covariate='p23650_i0',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step2_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='2')])

    # prepped-for time - step 3
    step3_var_explained_before, step3_residuals, step3_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step2_residuals,
                                                                                         covariate='prepped_for_time',
                                                                                         metabolite=metabolite,
                                                                                         covariate_type='continuous')

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step3_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='3')])

    # date at centre - step 4
    step4_var_explained_before, step4_residuals, step4_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step3_residuals,
                                                                                         covariate='p53_i0',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step4_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='4')])

    # date sample prepared - step 5
    step5_var_explained_before, step5_residuals, step5_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step4_residuals,
                                                                                         covariate='sample_prepared_date',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step5_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='5')])

    # date sample measured - step 6
    step6_var_explained_before, step6_residuals, step6_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step4_residuals,
                                                                                         covariate='sample_measured_date',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=step6_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='6')])

    # well - step 7
    step7_var_explained_before, final_residuals, step7_var_explained_after = regress_out(covariate_df=covariates,
                                                                                         data_df=step6_residuals,
                                                                                         covariate='p23660_i0',
                                                                                         metabolite=metabolite)

    vars_explained_by_all_covars_at_each_step = concat([vars_explained_by_all_covars_at_each_step,
                                                        var_explained_by_each_covar(covariate_df=covariates,
                                                                                    data_df=final_residuals,
                                                                                    covariates=covariate_info,
                                                                                    metabolite=metabolite,
                                                                                    stage='after (7)')])

    all_final_residuals = merge(all_final_residuals, final_residuals, on='eid', how='left')

    # bar.update()

all_final_residuals.iloc[:, 1:] = robust_scale(all_final_residuals.iloc[:, 1:])
all_final_residuals.to_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_post_regressions_robust_scaled_data.csv', index=False)
vars_explained_by_all_covars_at_each_step.to_csv('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/effects_during_regressions.csv', index=False)

# plot covariate effects over 'time'
vars_explained_by_all_covars_at_each_step['perc_var_explained'] = vars_explained_by_all_covars_at_each_step['var_explained'].apply(lambda x: x*100)
fig, ax = plt.subplots(figsize=[11.693, 8.268])
sns.lineplot(vars_explained_by_all_covars_at_each_step, x='stage', y='perc_var_explained', hue='covar_name', estimator='mean', errorbar='sd')
ax.set_ylabel('% variance explained')
ax.set_title('% variance explained by each covariate after each covariate is successively regressed out')
plt.legend(title='covariate')
fig.tight_layout()
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/effects_during_regressions.png')
fig.savefig('~/ch4_other_data_layers/4.1 metabolomics/4.1.1 covariates/plots/effects_during_regressions.svg')
