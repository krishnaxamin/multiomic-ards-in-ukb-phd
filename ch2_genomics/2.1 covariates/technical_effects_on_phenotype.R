' Getting technical covariate effects on phenotype, using logistic regression '

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

pan_ukbb_eur_eids <- read.csv('~/data/internal/cohort_eids/pan-ukbb-eur_eids.tsv', sep = '\t', header=FALSE)$V1
high_missingness_eids <- read.csv('~/data/internal/genomics/step4_high_missingness_samples_panukbb_eur.tsv', sep = '\t', header=FALSE)$V1

genomic_technical_covars <- read.csv('~/data/internal/genomics/genomics_covars.csv') %>%
  select(-p52, -p34) %>%
  filter(eid %in% pan_ukbb_eur_eids) %>%
  filter(!(eid %in% high_missingness_eids)) %>%
  mutate(plate_row = substr(p22008, start = 1, stop = 1)) %>%
  mutate(plate_column = substr(p22008, start = 2, stop= 3))

qu10_50_disease_fields <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')$disease_field
disease <- qu10_50_disease_fields[task_number]
phenotype <- read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv')[c('eid', disease)]


var_explained_by_covar <- function(covars_df, response_var_df, covar, response_var,
                                   covar_type = 'categorical', response_var_type = 'continuous') {
  
  covar_response_var = merge(covars_df[c('eid', covar)],
                             response_var_df, by = 'eid')
  
  colnames(covar_response_var) <- c('eid', 'covar', 'response_var')
  
  if (response_var_type == 'continuous'){
    if (covar_type == 'categorical'){
      model <- lm(response_var ~ factor(covar), data = covar_response_var)
    }
    else if (covar_type == 'continuous'){
      model <- lm(response_var ~ covar, data = covar_response_var)
    }
  }
  
  if (response_var_type == 'categorical'){
    if (covar_type == 'categorical'){
      model <- glm(response_var ~ factor(covar), data = covar_response_var, family = 'binomial')
    }
    else if (covar_type == 'continuous'){
      model <- glm(response_var ~ covar, data = covar_response_var, family = 'binomial')
    }
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
  
  # var_explained <- model_summary$r.squared
  
  return(list(coefficient = coefficient, pvalue = pvalue))
  
}


# Define the covariate_info list
covariate_info <- list(
  list(field = 'p22000', name = 'Batch', type = 'categorical'),
  list(field = 'p53_i0', name = 'Date', type = 'categorical'),
  list(field = 'p22008', name = 'Well', type = 'categorical'),
  list(field = 'p22007', name = 'Plate', type = 'categorical'),
  list(field = 'plate_row', name = 'Plate row', type = 'categorical'),
  list(field = 'plate_column', name = 'Plate column', type = 'categorical'),
  list(field = 'p31', name = 'Sex', type = 'categorical'),
  list(field = 'p21003_i0', name = 'Age', type = 'continuous'),
  list(field = 'p54_i0', name = 'Centre', type = 'categorical')
)
names(covariate_info) <- sapply(covariate_info, function(x) x$field)

phenotype_vars_explained <- data.frame()

for (x in covariate_info){
  print(x$name)  # flush these to output log output
  flush.console()
  Sys.sleep(1)
  start_time <- Sys.time()
  model_result <- var_explained_by_covar(covars_df = genomic_technical_covars,
                                         response_var_df = phenotype,
                                         covar = x$field,
                                         response_var = disease,
                                         covar_type = x$type,
                                         response_var_type = 'categorical')
  # if (is.null(model_result$var_explained)){
  #  model_result$var_explained <- -1
  # }
  model_result$covar <- x$name
  
  model_result_df <- data.frame(model_result)
  
  print(Sys.time() - start_time)
  flush.console()
  Sys.sleep(1)
  
  # write_csv(model_result_df, paste0('covariate_effects/', disease, '_', x$field, '_covariate_effects_on_phenotype.csv'))
  
  phenotype_vars_explained <- rbind(phenotype_vars_explained, model_result_df)
}

write_csv(phenotype_vars_explained, paste0('~/ch2_genomics/2.1 covariates/technical_effects/', disease, '_technical_effects_on_phenotype.csv'))
