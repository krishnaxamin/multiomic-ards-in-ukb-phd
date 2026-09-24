' Getting covariate effects on proteomics prior phenotypes, using logistic regression '

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

participant_covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv') %>%
  select(eid, p21003_i0, p31, p53_i0, p54_i0, p30901_i0, p30902_i0, Batch, plate_row, plate_column)

protein_covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv') %>%
  select(PlateID, Processing_StartDate)

protein_covars_mean_dates <- protein_covars %>%
  group_by(PlateID) %>%
  summarize(mean_processing_start_date = as.Date(mean(as.numeric(as.Date(Processing_StartDate))), origin = "1970-01-01")) %>%
  rename(p30901_i0 = PlateID)

protein_covars_dates_per_plate <- protein_covars %>%
  group_by(PlateID) %>%
  summarize(dates = paste0(unique(Processing_StartDate), collapse = '+')) %>%
  rename(p30901_i0 = PlateID)

covars <- merge(participant_covars, protein_covars_mean_dates, by = 'p30901_i0') %>%
  mutate(storage_time = as.Date(mean_processing_start_date) - as.Date(p53_i0))

qu10_50_disease_fields <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')$disease_field
disease <- qu10_50_disease_fields[task_number]
phenotype <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv')[c('eid', disease)]


var_explained_by_covar <- function(covars_df, response_var_df, covar, response_var,
                                   covar_type = 'categorical', response_var_type = 'continuous') {
  
  covar_response_var = merge(covars_df[c('eid', covar)],
                             response_var_df, by = 'eid')
  
  colnames(covar_response_var) <- c('eid', 'covar', 'response_var')
  
  # linear regression
  if (response_var_type == 'continuous'){
    if (covar_type == 'categorical'){
      model <- lm(response_var ~ factor(covar), data = covar_response_var)
    } else if (covar_type == 'continuous'){
      model <- lm(response_var ~ covar, data = covar_response_var)
    }
    
    model_summary <- summary(model)
    
    # multiple coeffs and pvalues returned for each category of categorical variable. return most significant results.
    coefficient_df <- data.frame(model_summary$coefficients[, 'Estimate'])[-c(1), , drop = FALSE]
    pvalue_df <- data.frame(model_summary$coefficients[, 'Pr(>|t|)'])[-c(1), , drop = FALSE]
    df <- cbind(coefficient_df, pvalue_df)
    colnames(df) <- c('log_odds', 'pvalue')
    df <- df %>% arrange(pvalue)
    coefficient <- df[1, 'log_odds']
    pvalue <- df[1, 'pvalue']
    var_explained <- model_summary$r.squared
  }
  
  # logistic regression
  if (response_var_type == 'categorical'){
    if (covar_type == 'categorical'){
      model <- glm(response_var ~ factor(covar), data = covar_response_var, family = 'binomial')
    } else if (covar_type == 'continuous'){
      model <- glm(response_var ~ covar, data = covar_response_var, family = 'binomial')
    }
    
    model_summary <- summary(model)
    
    # multiple coeffs and pvalues returned for each category of categorical variable. return most significant results.
    coefficient_df <- data.frame(model_summary$coefficients[, 'Estimate'])[-c(1), , drop = FALSE]
    pvalue_df <- data.frame(model_summary$coefficients[, 'Pr(>|z|)'])[-c(1), , drop = FALSE]
    df <- cbind(coefficient_df, pvalue_df)
    colnames(df) <- c('log_odds', 'pvalue')
    df <- df %>% arrange(pvalue)
    coefficient <- df[1, 'log_odds']
    pvalue <- df[1, 'pvalue']
    var_explained <- -1
    
  }
  
  
  return(list(coefficient = coefficient, pvalue = pvalue, var_explained = var_explained))
  
}


# Define the covariate_info list
covariate_info <- list(
  list(field = 'p21003_i0', name = 'age', type = 'continuous'),
  list(field = 'p31', name = 'sex', type = 'categorical'),
  list(field = 'p53_i0', name = 'date_centre', type = 'categorical'),
  list(field = 'p54_i0', name = 'centre', type = 'categorical'),
  list(field = 'p30901_i0', name = 'plate', type = 'categorical'),
  list(field = 'p30902_i0', name = 'well', type = 'categorical'),
  list(field = 'plate_row', name = 'plate_row', type = 'categorical'),
  list(field = 'plate_column', name = 'plate_column', type = 'categorical'),
  list(field = 'mean_processing_start_date', name = 'mean_date_processing', type = 'categorical'),
  list(field = 'storage_time', name = 'storage_time', type = 'continuous'),
  list(field = 'Batch', name = 'batch', type = 'categorical')
)
names(covariate_info) <- sapply(covariate_info, function(x) x$field)

phenotype_vars_explained <- data.frame()

if (var(phenotype[disease]) != 0){
  for (x in covariate_info){
    print(x$name)  # flush these to output log output
    flush.console()
    Sys.sleep(1)
    start_time <- Sys.time()
    model_result <- var_explained_by_covar(covars_df = covars,
                                           response_var_df = phenotype,
                                           covar = x$field,
                                           response_var = disease,
                                           covar_type = x$type,
                                           response_var_type = 'categorical')
  
    model_result$covar <- x$name
    
    model_result_df <- data.frame(model_result)
    
    print(Sys.time() - start_time)
    flush.console()
    Sys.sleep(1)
    
    # write_csv(model_result_df, paste0('covariate_effects/', disease, '_', x$field, '_covariate_effects_on_phenotype.csv'))
    
    phenotype_vars_explained <- rbind(phenotype_vars_explained, model_result_df)
  }
  
  write_csv(phenotype_vars_explained, paste0('~/ch3_proteomics/3.1 covariates/technical_effects/on_prior_phenotype/', disease, '_technical_effects_on_phenotype.csv'))
}
