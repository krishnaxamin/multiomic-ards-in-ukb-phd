""" Obtain per-disease age-of-onset profile and related data for the genomics cohort. """
from pandas import read_csv, DataFrame, concat, merge
from datetime import datetime
from collections import Counter
from math import ceil

import os


def common_unisex_disease_check(input_df, sex_df, cohort_num_ppl=1000, case_control_max=199,
                                cohort_num_males=10, cohort_num_females=10, unisex_max_rate=0.001, unisex_check=True):
    """
    Check that disease passes a certain case:control ratio and is common above a certain rate in both males and females.
    If unisex_check=False, e.g. when looking at diseases in males only, then the unisex check is bypassed.
    :param input_df:
    :param cohort_num_ppl:
    :param case_control_max:
    :param cohort_num_males:
    :param cohort_num_females:
    :param unisex_max_rate:
    :param unisex_check:
    :return:
    """

    input_df_with_sex = merge(input_df, sex_df, on='eid')

    disease_sex_counter = Counter(input_df_with_sex['sex'].tolist())
    if unisex_check:
        if disease_sex_counter['Male'] >= cohort_num_males * unisex_max_rate and disease_sex_counter['Female'] >= cohort_num_females * unisex_max_rate:
            disease_is_unisex = True
        else:
            disease_is_unisex = False
    else:
        disease_is_unisex = True

    disease_num_ppl = len(input_df_with_sex)
    if (cohort_num_ppl - disease_num_ppl) / disease_num_ppl <= case_control_max:
        disease_is_common = True
    else:
        disease_is_common = False

    disease_is_common_and_unisex = disease_is_unisex & disease_is_common

    return disease_is_common_and_unisex


def make_age_of_onset_time_series(input_df, max_age=100):
    """
    Base function for making time series (disease onset rate vs rounded age-of-onset).
    :param input_df:
    :param max_age:
    :return:
    """
    counter_dict = Counter(input_df['age_of_onset_rounded'].to_list())
    min_age = 0
    max_age_ceil = ceil(max_age)
    age_range = [i / 10.0 for i in range(min_age * 10, (max_age_ceil * 10), 5)]
    age_range.append(float(max_age_ceil))
    disease_onset_rate_list = [0.0] * len(age_range)
    num_diagnosed_list = [0.0] * len(age_range)
    for age in input_df['age_of_onset_rounded'].unique().tolist():
        disease_onset_rate_list[age_range.index(age)] = counter_dict[age] / len(input_df)
        num_diagnosed_list[age_range.index(age)] = counter_dict[age]

    disease_time_series_df = DataFrame({'age': age_range,
                                        'num_diagnosed': num_diagnosed_list,
                                        'disease_onset_rate': disease_onset_rate_list})
    return disease_time_series_df


def get_age_of_onset_info_for_common_unisex_diseases(subset_phenotypes, subset_eid_sex,
                                                     disease_codes_date_fields_list_func=None,
                                                     covariates_func=None,
                                                     case_control_ratio=199,
                                                     unisex_rate=0.001,
                                                     unisex_check_bool=True):
    """
    This inputs a subset population, which is represented by the phenotypes and sex information of that subset only.
    Then, it checks whether diseases are common and unisex within that population subset.
    Then, for each commmon-and-unisex disease, the age-of-onsets for each participant in this population subset is obtained.
    This returns a dictionary with elements:
    - 'age_of_onset': dictionary containing per-participant age-of-onsets for each disease {<disease>: <df>, ...}
    - 'max_age_of_onset': the maximum age_of_onset found for this subset of population across all diseases <int>
    :param subset_phenotypes:
    :param subset_eid_sex:
    :param disease_codes_date_fields_list_func:
    :param covariates_func:
    :param case_control_ratio:
    :param unisex_rate:
    :param unisex_check_bool:
    :return:
    """

    # check that same EIDs are contained within phenotype and sex information of the subset population
    assert(((subset_phenotypes.sort_values(by='eid')['eid'].reset_index(drop=True) == subset_eid_sex.sort_values(by='eid')['eid'].reset_index(drop=True)).sum() == len(subset_phenotypes)) &
           (len(subset_phenotypes) == len(subset_eid_sex)))

    disease_df_dict = {}  # dict populated with a df per disease. the df is keyed with the disease field. the df contains age-of-onset for each afflicted.
    max_age_of_onset_list = []
    onset_relative_to_assessment_df = DataFrame(columns=['disease_field',
                                                         'num_onset_before_assessment',
                                                         'num_onset_after_assessment'])

    subset_sexes_counts = Counter(subset_eid_sex['sex'].tolist())
    subset_num_males = subset_sexes_counts['Male']
    subset_num_females = subset_sexes_counts['Female']

    for disease_code_date_field in disease_codes_date_fields_list_func:
        # print(disease_code_date_field, flush=True)
        disease_code_source_field = 'p' + str(int(disease_code_date_field.split('p')[1]) + 1)
        disease_df = subset_phenotypes[
            subset_phenotypes[disease_code_date_field].notnull()][
            ['eid', disease_code_date_field, disease_code_source_field]].reset_index(drop=True)
        # checks if there's anyone with the disease
        if disease_df.empty:
            continue
        else:
            # checks if the disease appears with a certain rate in the cohort in both males and females, e.g. >1 in every 1000 males and >1 in every 1000 females
            if not common_unisex_disease_check(disease_df,
                                               sex_df=subset_eid_sex,
                                               cohort_num_ppl=len(subset_phenotypes),
                                               case_control_max=case_control_ratio,
                                               cohort_num_males=subset_num_males,
                                               cohort_num_females=subset_num_females,
                                               unisex_max_rate=unisex_rate,
                                               unisex_check=unisex_check_bool):
                continue

        # add p52, p34, p53_i0 to disease_df
        disease_df = merge(disease_df, covariates_func[['eid', 'p52', 'p34', 'p53_i0']], on='eid')

        age_of_onset_list = []
        age_of_onset_rounded_list = []
        onset_relative_to_assessment_list = []
        for i in range(len(disease_df)):

            # make D.O.B from p52, p34 - this could be done also on the OG imported df and replace p52 and p34
            # accounts for data-coding 819 (dates that have other meanings)
            # (see UKB data-coding 819 - these pseudodates represent errors/oddities in the first occurrence dates)
            if any(x in disease_df[disease_code_date_field][i] for x in ['1900-01-01', '1901-01-01', '2037-07-07',
                                                                         'Code has no event date',
                                                                         'Code has event date before participant\'s date of birth',
                                                                         'Code has event date in the future and is presumed to be a place-holder or other system default']):
                disease_df = disease_df.drop(i)
                continue
            # if disease_df[disease_code_date_field][i] in ['1900-01-01', '1901-01-01', '2037-07-07',
            # 'Code has no event date', 'Code has event date before participant\'s date of birth',
            # 'Code has event date in the future and is presumed to be a place-holder or other system default']:
            #     disease_df = disease_df.drop(i)
            #     continue
            elif any(x in disease_df[disease_code_date_field][i] for x in ['1902-02-02', '1903-03-03',
                                                                           'Code has event date matching participant\'s date of birth',
                                                                           'Code has event date after participant\'s date of birth and falls in the same calendar year as date of birth']):
                age_of_onset_list.append(0.0)
                age_of_onset_rounded_list.append(0.0)
                onset_relative_to_assessment_list.append('before')

            else:
                dob = datetime.strptime('1 ' + disease_df['p52'][i] + ' ' + str(disease_df['p34'][i]), '%d %B %Y')
                date_of_assessment = datetime.strptime(disease_df['p53_i0'][i], '%Y-%m-%d')
                date_of_onset = datetime.strptime(disease_df[disease_code_date_field][i], '%Y-%m-%d')

                # get age-of-onset from D.O.B and date-of-onset
                if date_of_onset.month < dob.month:
                    age_of_onset = (date_of_onset.year - dob.year - 1) + (12 - (dob.month - date_of_onset.month)) / 12
                else:
                    age_of_onset = (date_of_onset.year - dob.year) + (date_of_onset.month - dob.month) / 12
                age_of_onset_list.append(age_of_onset)

                if age_of_onset % 1 == 0.25:
                    age_of_onset_rounded = age_of_onset + 0.25
                else:
                    age_of_onset_rounded = round(age_of_onset * 2) / 2
                age_of_onset_rounded_list.append(age_of_onset_rounded)

                # onset before or after assessment
                if date_of_onset > date_of_assessment:
                    onset_relative_to_assessment = 'after'
                elif date_of_onset == date_of_assessment:
                    onset_relative_to_assessment = 'on'
                else:
                    onset_relative_to_assessment = 'before'
                onset_relative_to_assessment_list.append(onset_relative_to_assessment)

        disease_df['age_of_onset'] = age_of_onset_list
        disease_df['age_of_onset_rounded'] = age_of_onset_rounded_list
        onset_relative_to_assessment_counter = Counter(onset_relative_to_assessment_list)
        onset_relative_to_assessment_df = concat(
            [onset_relative_to_assessment_df, DataFrame([{'disease_field': disease_code_date_field,
                                                          'num_onset_before_assessment':
                                                              onset_relative_to_assessment_counter['before'],
                                                          'num_onset_after_assessment':
                                                              onset_relative_to_assessment_counter['after']}])])

        disease_df = disease_df[['eid', 'age_of_onset', 'age_of_onset_rounded']]
        disease_df_dict[disease_code_date_field] = disease_df
        max_age_of_onset_list.append(max(disease_df['age_of_onset_rounded']))

    return {'age_of_onsets': disease_df_dict,
            'onset_relative_to_assessment': onset_relative_to_assessment_df,
            'max_age_of_onset': max(max_age_of_onset_list)}


""" Load data """
phenotypes = read_csv('~/data/external/genomics/array-genotyping_phenotype_file.csv')
covariates = read_csv('~/data/external/genomics/genomics_covars.csv')
pan_ukbb_eur_eids = list(read_csv('~/data/internal/cohort_eids/genomics_pan-ukbb-eur_eids.tsv', delim_whitespace=True, header=None)[0])
high_missingness_eids = list(read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
                                      delim_whitespace=True, header=None)[0])

filtered_phenotypes = phenotypes[(phenotypes['eid'].isin(pan_ukbb_eur_eids)) & (~phenotypes['eid'].isin(high_missingness_eids))]
filtered_covariates = covariates[(covariates['eid'].isin(pan_ukbb_eur_eids)) & (~covariates['eid'].isin(high_missingness_eids))]

eid_sex = filtered_covariates[['eid', 'p31']]
eid_sex.columns = ['eid', 'sex']

""" Get age-of-onset for each person with a disease """

# get fields which contain first occurrence dates for each disease code
disease_codes_date_fields_list = [x for x in filtered_phenotypes.columns.to_list()[1:] if int(x.split('_')[0].split('p')[1]) % 2 == 0 and 'p13' in x]

per_disease_age_of_onset_info = {}
per_disease_onset_relative_to_assessment_info = {}
max_age_of_onset_list = []

# everyone
print('Age-of-onset: everyone', flush=True)
everyone_info = get_age_of_onset_info_for_common_unisex_diseases(subset_phenotypes=filtered_phenotypes,
                                                                 subset_eid_sex=eid_sex,
                                                                 disease_codes_date_fields_list_func=disease_codes_date_fields_list,
                                                                 covariates_func=filtered_covariates)
per_disease_age_of_onset_info['everyone'] = everyone_info['age_of_onsets']
per_disease_onset_relative_to_assessment_info['everyone'] = everyone_info['onset_relative_to_assessment']
max_age_of_onset_list.append(everyone_info['max_age_of_onset'])

# by sex - unisex_check_bool=False since only looking at one sex
sexes = list(set(eid_sex.sex))
for sex in sexes:
    print('Age-of-onset: ' + sex, flush=True)
    eid_sex_one_sex = eid_sex[eid_sex['sex'] == sex].copy()
    eids_one_sex = list(eid_sex_one_sex['eid'])
    phenotypes_one_sex = filtered_phenotypes[filtered_phenotypes['eid'].isin(eids_one_sex)].copy()
    info_one_sex = get_age_of_onset_info_for_common_unisex_diseases(subset_phenotypes=phenotypes_one_sex,
                                                                    subset_eid_sex=eid_sex_one_sex,
                                                                    disease_codes_date_fields_list_func=disease_codes_date_fields_list,
                                                                    covariates_func=filtered_covariates,
                                                                    unisex_check_bool=False)
    per_disease_age_of_onset_info[sex] = info_one_sex['age_of_onsets']
    per_disease_onset_relative_to_assessment_info[sex] = info_one_sex['onset_relative_to_assessment']
    max_age_of_onset_list.append(info_one_sex['max_age_of_onset'])

""" Export common and unisex disease fields in the genomics cohort (everyone, by sex, by ancestry group) """
common_unisex_diseases = DataFrame({'disease_field': disease_codes_date_fields_list})
print('Exporting common_unisex_diseases', flush=True)
for key in per_disease_age_of_onset_info.keys():
    specific_common_unisex_diseases = list(per_disease_age_of_onset_info[key].keys())
    specific_common_unisex_disease_df = DataFrame({'disease_field': specific_common_unisex_diseases,
                                                   key.lower(): [1] * len(specific_common_unisex_diseases)})
    common_unisex_diseases = merge(common_unisex_diseases, specific_common_unisex_disease_df, on='disease_field', how='outer')

common_unisex_diseases = common_unisex_diseases[common_unisex_diseases.iloc[:, 1:].sum(axis=1) > 0]
common_unisex_diseases.to_csv('~/data/internal/genomics/genomics_qc2_pan_ukbb_eur_common_unisex_disease_date_fields.csv', index=False)

""" Get maximum age-of-onset across all cases - this sets the upper bound of 'age' for easy plotting """
max_age_of_onset = max(max_age_of_onset_list)

""" Create disease-specific age-of-onset profiles """
for key in per_disease_age_of_onset_info.keys():
    print('Creating disease-specific age-of-onset profiles: ' + key, flush=True)
    specific_disease_df_dict = per_disease_age_of_onset_info[key]
    for disease in specific_disease_df_dict.keys():
        time_series = make_age_of_onset_time_series(specific_disease_df_dict[disease],
                                                    max_age=max_age_of_onset)
        export_path_dir = '~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_' + key.lower()
        if not os.path.exists(export_path_dir):
            os.makedirs(export_path_dir)
        time_series.to_csv(export_path_dir + '/' + disease + '_' + key.lower() + '_time_series.csv', index=False)
