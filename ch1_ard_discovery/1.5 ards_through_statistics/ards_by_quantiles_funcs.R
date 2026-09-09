library(tidyverse)

' Functions used to extract age-of-onset quantiles and use those to identify 
age-related diseases. '

# this function currently only valid for genomics cohort - likely to be deprecated if data for genomics cohort is reworked to be like those for the other omics cohorts
age_of_onset_stats_extractor_no_num_diagnosed_in_time_series <- function(quantile_probs){
  
  time_series_file_names <- as.list(list.files(path = '~/data/internal/time_series/genomics_qcv1_time_series/raw/', pattern = '*_time_series.csv'))
  num_time_series <- length(time_series_file_names)
  disease_field_chapters <- read.csv('~/data/internal/phenotype_coding/disease_fields_icd10_chapters.csv')
  disease_fields_to_icd10 <- read.csv('~/data/internal/phenotype_coding/disease_fields_icd10_info.csv')
  
  num_ppl_per_disease_df <- read.csv('~/data/internal/genomics/genomics_qc1_common-unisex-disease_num-diagnosed_info-content_case-control.csv')

  stats_df <- data.frame(disease = character(),
                         quantile = numeric())
  
  # Get age-of-onset stats for each disease
  for(i in seq_along(l=time_series_file_names)) {
    time_series <- read.csv(paste('~/data/internal/time_series/genomics_qcv1_time_series/raw/', time_series_file_names[[i]], sep=''))
    disease_field <- str_split(time_series_file_names[[i]], '_')[[1]][1]
    icd10_chapter <- disease_field_chapters[disease_field_chapters$disease_field == disease_field, "icd10_chapter"]
    
    num_cases <- num_ppl_per_disease_df[num_ppl_per_disease_df$disease == disease_field, 'num_cases']
    time_series$onset_counts <- round(time_series$disease_onset_rate * num_cases)
    onset_counts_spelt_out <- rep(time_series$age, time_series$onset_counts)
    disease_quantile <- quantile(onset_counts_spelt_out, probs = quantile_probs)[[1]]
    
    disease_age_of_onset_df <- data.frame(disease_field = disease_field,
                                          quantile = disease_quantile)
    
    stats_df <- rbind(stats_df, disease_age_of_onset_df)
  }
  
  merge(merge(stats_df, disease_fields_to_icd10, by = 'disease_field'), num_ppl_per_disease_df, by = 'disease_field')
  
}

# this function valid for time series with num_diagnosed field. specificity = 'everyone', 'male', 'female'
age_of_onset_stats_extractor_num_diagnosed_in_time_series <- function(quantile_probs, omic, specificity){
  
  path <- paste0(omic, '_pan_ukbb_eur_time_series/common_unisex_', tolower(specificity))
  
  time_series_file_names <- as.list(list.files(path = paste0('~/data/internal/time_series/', path, '/'), pattern = '*_time_series.csv'))
  num_time_series <- length(time_series_file_names)
  disease_field_chapters <- read.csv('~/data/internal/phenotype_coding/disease_fields_icd10_chapters.csv')
  disease_fields_to_icd10 <- read.csv('~/data/internal/phenotype_coding/disease_fields_icd10_info.csv')

  stats_df <- data.frame(disease = character(),
                         quantile = numeric())
  
  # Get age-of-onset stats for each disease
  for(i in seq_along(l=time_series_file_names)) {
    time_series <- read.csv(paste0('~/data/internal/time_series/', path, '/', time_series_file_names[[i]]))
    disease_field <- str_split(time_series_file_names[[i]], '_')[[1]][1]
    icd10_chapter <- disease_field_chapters[disease_field_chapters$disease_field == disease_field, "icd10_chapter"]
    
    onset_counts_spelt_out <- rep(time_series$age, time_series$num_diagnosed)
    disease_quantile <- quantile(onset_counts_spelt_out, probs = quantile_probs)[[1]]
    
    disease_age_of_onset_df <- data.frame(disease_field = disease_field,
                                          quantile = disease_quantile)
    
    stats_df <- rbind(stats_df, disease_age_of_onset_df)
  }
  
  merge(merge(stats_df, disease_fields_to_icd10, by = 'disease_field'), 
        read.csv(paste0('~/data/internal/', omic, '/', omic, '_pan_ukbb_eur_common-unisex-disease_num-diagnosed_case-control.csv')), by = 'disease_field')
  
}
