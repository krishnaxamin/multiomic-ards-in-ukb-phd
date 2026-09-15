"""
Identifies diseases that occur in more than 1/1000 males AND 1/1000 females (unisex).
Builds time series of disease onset rate vs age for each identified disease.
Inputs cohort browser data downloaded from the UKB on DNAnexus using command in cohort_data_download.txt.
"""

from pandas import read_csv, DataFrame, concat
from datetime import datetime
from collections import Counter
from math import ceil


# check that disease is common above a set rate in both males and females
def common_unisex_disease_check(input_df, cohort_num_ppl=1000, case_control_max=199,
                                cohort_num_males=10, cohort_num_females=10, unisex_max_rate=0.001):
    disease_sex_counter = Counter(input_df['p31'].tolist())
    if disease_sex_counter['Male'] >= cohort_num_males * unisex_max_rate and disease_sex_counter['Female'] >= cohort_num_females * unisex_max_rate:
        disease_is_unisex = True
    else:
        disease_is_unisex = False

    disease_num_ppl = len(input_df)
    if (cohort_num_ppl - disease_num_ppl) / disease_num_ppl <= case_control_max:
        disease_is_common = True
    else:
        disease_is_common = False

    disease_is_common_and_unisex = disease_is_unisex & disease_is_common

    return disease_is_common_and_unisex


# make time series (disease onset rate vs rounded age-of-onset)
def make_age_of_onset_time_series(input_df, max_age=100):
    counter_dict = Counter(input_df['age_of_onset_rounded'].to_list())
    min_age = 0
    max_age_ceil = ceil(max_age)
    age_range = [i / 10.0 for i in range(min_age * 10, (max_age_ceil * 10), 5)]
    age_range.append(float(max_age_ceil))
    disease_onset_rate_list = [0.0] * len(age_range)
    for age in input_df['age_of_onset_rounded'].unique().tolist():
        disease_onset_rate_list[age_range.index(age)] = counter_dict[age] / len(input_df)

    disease_time_series_df = DataFrame({'age': age_range, 'disease_onset_rate': disease_onset_rate_list})
    return disease_time_series_df


# makes times series, for each sex
def make_age_of_onset_time_series_by_sex(input_df, max_age=100):
    by_sex_counter = Counter(input_df['sex'].to_list())
    disease_time_series_by_sex_df_list = []
    for sex in list(by_sex_counter.keys()):
        input_df_sex = input_df[input_df['sex'] == sex]
        disease_time_series_df = make_age_of_onset_time_series(input_df_sex, max_age=max_age)
        disease_time_series_df['sex'] = sex
        disease_time_series_by_sex_df_list.append(disease_time_series_df)
    return concat(disease_time_series_by_sex_df_list)


df = read_csv('~/data/external/genomics/array-genotyping_phenotype_file.csv')

# number of males/females in the cohort
sex_counter = Counter(df['p31'].tolist())
num_males = sex_counter['Male']
num_females = sex_counter['Female']
num_ppl = len(df)

# get fields which contain first occurrence dates for each disease code
disease_codes_date_fields_list = [x for x in df.columns.to_list()[1:] if int(x.split('_')[0].split('p')[1]) % 2 == 0 and 'p13' in x]

# see UKB data-coding 819 - these pseudodates represent errors/oddities in the first occurrence dates
# the dates are now directly referenced when they are needed rather than being stored in variables which are referenced only once
# data_coding_819_exclusions = ['1900-01-01', '1901-01-01', '2037-07-07']
# data_coding_819_age_of_onset_is_zero = ['1902-02-02', '1903-03-03']

# loop returns a dictionary of dataframes, with one df per disease and dfs containing age-of-onset info
disease_df_dict = {}
disease_time_series_dict = {}
max_age_of_onset_list = []
onset_relative_to_assessment_df = DataFrame(columns=['disease_code',
                                                     'num_onset_before_assessment',
                                                     'num_onset_after_assessment'])
for disease_code_date_field in disease_codes_date_fields_list:
    # print(disease_code_date_field)
    # disease_code_date_field = 'p130016'
    disease_code_source_field = 'p' + str(int(disease_code_date_field.split('p')[1]) + 1)
    disease_df = df[df[disease_code_date_field].notnull()][['eid', 'p31', 'p21000_i0', 'p52', 'p34', 'p53_i0',
                                                          disease_code_date_field,
                                                          disease_code_source_field]].reset_index(drop=True)
    # checks if there's anyone with the disease
    if disease_df.empty:
        continue
    else:
        # checks if the disease appears with a certain rate in the cohort in both males and females, e.g. >1 in every 1000 males and >1 in every 1000 females
        if not common_unisex_disease_check(disease_df, cohort_num_ppl=num_ppl, case_control_max=199,
                                           cohort_num_males=num_males, cohort_num_females=num_females,
                                           unisex_max_rate=0.001):
            continue

    age_of_onset_list = []
    age_of_onset_rounded_list = []
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
        elif any(x in disease_df[disease_code_date_field][i] for x in ['1902-02-02', '1903-03-03',
                                                                       'Code has event date matching participant\'s date of birth',
                                                                       'Code has event date after participant\'s date of birth and falls in the same calendar year as date of birth']):
            age_of_onset_list.append(0.0)
            age_of_onset_rounded_list.append(0.0)
        else:
            dob = datetime.strptime('1 ' + disease_df['p52'][i] + ' ' + str(disease_df['p34'][i]), '%d %B %Y')
            date_of_assessment = datetime.strptime(disease_df['p53_i0'][i], '%Y-%m-%d')
            date_of_onset = datetime.strptime(disease_df[disease_code_date_field][i], '%Y-%m-%d')

            # get age-of-onset from D.O.B and date-of-onset
            if date_of_onset.month < dob.month:
                age_of_onset = (date_of_onset.year - dob.year - 1) + (12 - (dob.month - date_of_onset.month))/12
            else:
                age_of_onset = (date_of_onset.year - dob.year) + (date_of_onset.month - dob.month)/12
            age_of_onset_list.append(age_of_onset)

            if age_of_onset % 1 == 0.25:
                age_of_onset_rounded = age_of_onset + 0.25
            else:
                age_of_onset_rounded = round(age_of_onset * 2)/2
            age_of_onset_rounded_list.append(age_of_onset_rounded)

    disease_df['age_of_onset'] = age_of_onset_list
    disease_df['age_of_onset_rounded'] = age_of_onset_rounded_list

    # disease_df['onset_relative_to_assessment'] = onset_relative_to_assessment_list
    # disease_df = disease_df[['eid', 'age_of_onset', 'age_of_onset_rounded', 'onset_relative_to_assessment']]
    disease_df = disease_df[['eid', 'p31', 'p21000_i0', 'age_of_onset', 'age_of_onset_rounded']]
    disease_df.columns = ['eid', 'sex', 'self_ethnicity', 'age_of_onset', 'age_of_onset_rounded']
    disease_df_dict[disease_code_date_field] = disease_df
    max_age_of_onset_list.append(max(disease_df['age_of_onset_rounded']))

# get the date fields for the selected diseases
common_unisex_disease_date_fields = DataFrame({'disease_field': list(disease_df_dict.keys())})
common_unisex_disease_date_fields.to_csv('~/data/internal/genomics/genomics_qc1_common_unisex_disease_date_fields.csv', index=False)

# gets max age-of-onset across all diseases
max_age_of_onset = max(max_age_of_onset_list)

# makes time series for each disease (raw; by sex; by self-reported ethnicity; by ancestry group (determined by
# self-reported ethnicity)) and exports it to a csv
for disease_taken_forward in list(disease_df_dict.keys()):
    print(disease_taken_forward)
    time_series = make_age_of_onset_time_series(disease_df_dict[disease_taken_forward],
                                                max_age=max_age_of_onset)
    disease_time_series_dict[disease_taken_forward] = time_series
    time_series.to_csv('~/data/internal/time_series/genomics_qcv1_time_series/raw/' + disease_taken_forward + '_time_series.csv',
                       index=False)

    time_series_by_sex = make_age_of_onset_time_series_by_sex(disease_df_dict[disease_taken_forward],
                                                              max_age=max_age_of_onset)
    time_series_by_sex.to_csv('~/data/internal/time_series/genomics_qcv1_time_series/by-sex/' + disease_taken_forward + '_by-sex_time_series.csv',
                              index=False)
