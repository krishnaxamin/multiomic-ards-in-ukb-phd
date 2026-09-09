"""
Rework the file mapping the UKB disease fields to ICD-10 chapters to include ICD-10 three letter codes and ICD-10
descriptions. The 'expanded' file includes sub-descriptions within an ICD-10 code,
e.g. different classifications of Cholera
"""

import pandas as pd

icd10_coding = pd.read_csv('~/data/external/phenotype_coding/icd10_codes-only_reworked_ukb-info.csv')
disease_field_to_icd10_chapter = pd.read_csv('~/data/external/phenotype_coding/disease_fields_icd10_chapters.csv')

icd10_chapters_codes_ukb = icd10_coding[['chapter', 'three_letter_code', 'in_ukb_first_occurrences']]

# code is appended twice to cover the two UKB fields that correspond to each disease (one field for date; one for source)
three_letter_code_list = []
for i in range(len(icd10_chapters_codes_ukb)):
    three_letter_code_list.append(icd10_chapters_codes_ukb['three_letter_code'][i])
    three_letter_code_list.append(icd10_chapters_codes_ukb['three_letter_code'][i])

disease_field_to_icd10_three_letter = disease_field_to_icd10_chapter
disease_field_to_icd10_three_letter['icd10_three_letter'] = three_letter_code_list

# maps UKB disease fields to ICD10 three-letter codes, with the matching ICD10 description
disease_field_to_icd10_three_letter_with_meaning_by_icd10_coding = disease_field_to_icd10_three_letter.merge(icd10_coding[['coding', 'meaning']], left_on='icd10_three_letter', right_on='coding').drop('coding', axis=1)

disease_field_to_icd10_three_letter_with_meaning_by_icd10_coding.to_csv('~/data/internal/phenotype_coding/disease_fields_icd10_info.csv', index=False)
