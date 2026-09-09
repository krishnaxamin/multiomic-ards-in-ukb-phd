' Getting lifestyle effects on genomics phenotype, using logistic regression '

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

covars <- read_csv('~/data/internal/genomics/genomics_lifestyles.csv', num_threads = 1, show_col_types = FALSE)

qu10_50_disease_fields <- read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt', num_threads = 1, show_col_types = FALSE)$disease_field
disease <- qu10_50_disease_fields[task_number]
phenotype <- read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv', num_threads = 1, show_col_types = FALSE)[c('eid', disease)]

# limit cohort to the chosen 425,223 people
pan_ukbb_eur_eids <- read_tsv('~/data/internal/cohort_eids/pan-ukbb-eur_eids.tsv', col_names = c('eid', 'eid1', 'x'),
                              num_threads = 1, show_col_types = FALSE)
step4_high_missigness_eids <- read_tsv('~/data/internal/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
                                       col_names = c('eid', 'eid1', 'x'),
                                       num_threads = 1, show_col_types = FALSE)
# keep pan-ukbb-eur eids, remove step-4 high-missingness eids
pan_ukbb_eur_covars <- covars %>%
  filter(eid %in% pan_ukbb_eur_eids$eid) %>%
  filter(!eid %in% step4_high_missigness_eids$eid)
pan_ukbb_eur_phenotype <- phenotype %>%
  filter(eid %in% pan_ukbb_eur_eids$eid) %>%
  filter(!eid %in% step4_high_missigness_eids$eid)

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


# Define the lifestyle_info list
lifestyle_info <- list(
  list(field = 'p1558_i0', name = 'Alcohol intake frequency', type = 'categorical'),
  list(field = 'p20116_i0', name = 'Smoking status', type = 'categorical'),
  list(field = 'p23099_i0', name = '% body fat', type = 'continuous'),
  list(field = 'p1160_i0', name = 'Hours slept per day', type = 'continuous'),
  list(field = 'p189', name = 'TDI', type = 'continuous'),
  list(field = 'p6138_i0', name = 'Education', type = 'categorical'),
  list(field = 'p884_i0', name = 'Days/week moderate activity', type = 'continuous'),
  list(field = 'grip_strength', name = 'Grip strength', type = 'continuous'),
  list(field = 'p1289_i0', name = 'Cooked veg intake', type = 'continuous'),
  list(field = 'p1299_i0', name = 'Raw veg intake', type = 'continuous'),
  list(field = 'p1309_i0', name = 'Fresh fruit intake', type = 'continuous'),
  list(field = 'p1319_i0', name = 'Dried fruit intake', type = 'continuous'),
  list(field = 'p1329_i0', name = 'Oily fish intake', type = 'categorical'),
  list(field = 'p1339_i0', name = 'Non-oily fish intake', type = 'categorical'),
  list(field = 'p1349_i0', name = 'Processed meat intake', type = 'categorical'),
  list(field = 'p1359_i0', name = 'Poultry intake', type = 'categorical'),
  list(field = 'p1369_i0', name = 'Beef intake', type = 'categorical'),
  list(field = 'p1379_i0', name = 'Lamb intake', type = 'categorical'),
  list(field = 'p1389_i0', name = 'Pork intake', type = 'categorical'),
  list(field = 'p1438_i0', name = 'Bread intake', type = 'continuous'),
  list(field = 'p1458_i0', name = 'Cereal intake', type = 'continuous'),
  list(field = 'p1528_i0', name = 'Water intake', type = 'continuous')
)
names(lifestyle_info) <- sapply(lifestyle_info, function(x) x$field)

phenotype_vars_explained <- data.frame()

for (x in lifestyle_info){
  print(x$name)  # flush these to output log output
  flush.console()
  Sys.sleep(1)
  start_time <- Sys.time()
  model_result <- var_explained_by_covar(covars_df = pan_ukbb_eur_covars,
                                         response_var_df = pan_ukbb_eur_phenotype,
                                         covar = x$field,
                                         response_var = disease,
                                         covar_type = x$type,
                                         response_var_type = 'categorical')
  
  model_result$covar <- x$name
  
  model_result_df <- data.frame(model_result)
  
  print(Sys.time() - start_time)
  flush.console()
  Sys.sleep(1)
  
  # write_csv(model_result_df, paste0('lifestyle_effects/', disease, '_', x$field, '_lifestyle_effects_on_phenotype.csv'))
  
  phenotype_vars_explained <- rbind(phenotype_vars_explained, model_result_df)
}

write_csv(phenotype_vars_explained, paste0('~/ch2_genomics/2.1 covariates/lifestyle_effects/', disease, '_lifestyle_effects_on_phenotype.csv'), num_threads = 1)
