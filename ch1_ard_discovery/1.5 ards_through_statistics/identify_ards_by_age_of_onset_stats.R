library(tidyverse)

source('~/ch1_ard_discovery/1.5 ards_through_statistics/ards_by_quantiles_funcs.R')

" Identify ARDs using statistics calculated on diseases' age-of-onset profiles, e.g. median > 50; 5% quantile > 50 "

# ARDs as qu.10 > 50 ------

## Genomics -------

# Genomics QC 1
genomics_qu10_stats <- age_of_onset_stats_extractor_no_num_diagnosed_in_time_series(0.1)
genomics_qu10_50 <- genomics_qu10_stats[genomics_qu10_stats$quantile > 50, ]
write_csv(genomics_qu10_50, '~/data/internal/phenotype_coding/genomics_qu10_50_icd10-info_num-diagnosed-info.csv')

# Genomics QC 2
genomics_qu10_stats <- age_of_onset_stats_extractor_num_diagnosed_in_time_series(0.1, omic = 'genomics', specificity = 'everyone')
genomics_qu10_50 <- genomics_qu10_stats[genomics_qu10_stats$quantile > 50, ]
write_csv(genomics_qu10_50, '~/data/internal/phenotype_coding/genomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')

## Metabolomics ------

metabolomics_qu10_stats <- age_of_onset_stats_extractor_num_diagnosed_in_time_series(0.1, omic = 'metabolomics', specificity = 'everyone')
metabolomics_qu10_50 <- metabolomics_qu10_stats[metabolomics_qu10_stats$quantile > 50, ]
write_csv(metabolomics_qu10_50, '~/data/internal/phenotype_coding/metabolomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')

## Proteomics

proteomics_qu10_stats <- age_of_onset_stats_extractor_num_diagnosed_in_time_series(0.1, omic = 'proteomics', specificity = 'everyone')
proteomics_qu10_50 <- proteomics_qu10_stats[proteomics_qu10_stats$quantile > 50, ]
write_csv(proteomics_qu10_50, '~/data/internal/phenotype_coding/proteomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')

## Comparing ------

genomics_qc1_qu10_50 <- read.csv('~/data/internal/phenotype_coding/genomics_qu10_50_icd10-info_num-diagnosed-info.csv')
metabolomics_qu10_50 <- read.csv('~/data/internal/phenotype_coding/metabolomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')
proteomics_qu10_50 <- read.csv('~/data/internal/phenotype_coding/proteomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')
genomics_qu10_50 <- read.csv('~/data/internal/phenotype_coding/genomics_pan_ukbb_eur_qu10_50_icd10-info_num-diagnosed-info.csv')

### Genomics vs genomics ------
# qu10_50 from genomics QC1 and QC2 share 69 ARDs
genomics_qc1_qc2_qu10_50 <- merge(genomics_qu10_50, genomics_qc1_qu10_50 %>%
                                         select(!info_content_nat),
                                       by = c('disease_field', 'category', 'icd10_chapter', 'icd10_three_letter', 'meaning'),
                                       suffixes = c('_genomics_qc2', '_genomics_qc1'))
for (x in genomics_qu10_50$disease_field){
  if (x %in% genomics_qc1_qc2_qu10_50$disease_field){
  }
  else{
    print(x)  # the extra one from QC2 is p131972
  }
}

### Genomics vs metabolomics ------
# 69 diseases: qu10_50 from metabolomics includes those from genomics QC1 qu10_50
metabolomics_genomics_qu10_50 <- merge(metabolomics_qu10_50, genomics_qu10_50 %>%
                                         select(!info_content_nat),
                                       by = c('disease_field', 'category', 'icd10_chapter', 'icd10_three_letter', 'meaning'),
                                       suffixes = c('_metabolomics', '_genomics'))

# 70 diseases: qu10_50 from metabolomics includes those from genomics QC2 qu10_50
metabolomics_genomics_qu10_50 <- merge(metabolomics_qu10_50, genomics_qu10_50,
                                       by = c('disease_field', 'category', 'icd10_chapter', 'icd10_three_letter', 'meaning'),
                                       suffixes = c('_metabolomics', '_genomics'))

### All omics ------
# 68 diseases: qu10_50 from proteomics, metabolomics and genomics QC2
all_omics_qu10_50 <- merge(proteomics_qu10_50, metabolomics_genomics_qu10_50,
                           by = c('disease_field', 'category', 'icd10_chapter', 'icd10_three_letter', 'meaning'))
all_omics_qu10_50 <- all_omics_qu10_50 %>% rename(quantile_proteomics = quantile,
                                                  num_cases_proteomics = num_cases,
                                                  case_control_ratio_proteomics = case_control_ratio)
write_csv(all_omics_qu10_50, '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')

### New ARDs vs old (from genomics QC1) ARDs ------
for (x in genomics_qc1_qu10_50$disease_field){
  if (x %in% all_omics_qu10_50$disease_field){
  }
  else{
    print(x)  # new ARDs lack p130660 and p132022 from old ARDs (made from genomics QC1)
  }
}

for (x in all_omics_qu10_50$disease_field){
  if (x %in% genomics_qc1_qu10_50$disease_field){
  }
  else{
    print(x)  # new ARDs gain p131972 from old ARDs
  }
}
