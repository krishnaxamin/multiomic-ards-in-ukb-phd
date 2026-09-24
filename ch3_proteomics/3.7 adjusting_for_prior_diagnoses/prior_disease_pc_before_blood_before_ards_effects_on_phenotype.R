' Getting previous disease PCx effects on proteomics phenotype, using logistic regression '

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

qu10_50_disease_fields <- read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt', num_threads = 1, show_col_types = FALSE)$disease_field
disease <- qu10_50_disease_fields[task_number]
phenotype <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv'), num_threads = 1, show_col_types = FALSE)[c('eid', disease)]
pcs <- read_csv(paste0('~/data/internal/proteomics/prior_disease_info/principal_components_before_blood_before_ards/', disease, '.csv'), num_threads = 1, show_col_types = FALSE)

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


# Define the pcs_info list
pcs_info <- lapply(1:(ncol(pcs) - 1), function(x) list(field = paste0('PC', x), name = paste0('PC', x), type = 'continuous'))
names(pcs_info) <- sapply(pcs_info, function(x) x$field)

phenotype_vars_explained <- data.frame()

if (var(phenotype[disease]) != 0){
  for (x in pcs_info){
    print(x$name)  # flush these to output log output
    flush.console()
    Sys.sleep(1)
    start_time <- Sys.time()
    model_result <- var_explained_by_covar(covars_df = pcs,
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
    
    # write_csv(model_result_df, paste0('lifestyle_effects/', disease, '_', x$field, '_lifestyle_effects_on_phenotype.csv'))
    
    phenotype_vars_explained <- rbind(phenotype_vars_explained, model_result_df)
  }
  
  write_csv(phenotype_vars_explained, paste0('ch3_proteomics/3.7 adjusting_for_prior_diagnoses/before_blood_before_ards_effect_on_phenotype/', disease, '_prior_disease_pc_effects_on_phenotype.csv'))
}
