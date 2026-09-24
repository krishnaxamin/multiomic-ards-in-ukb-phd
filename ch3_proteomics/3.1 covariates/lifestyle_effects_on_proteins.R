' Getting lifestyle effects on proteomics data, using linear regression.
Script enabled for use on pre- and post-regressions data.'

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])


# load in specific protein data
proteins <- read.csv('~/data/internal/proteomics/protein_fields.csv')$protein_field
protein <- proteins[task_number]
protein_data <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_data.csv')[c('eid', protein)]

covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv')

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
    # beta not required as only looking at R^2
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
  
  
  return(list(coefficient = coefficient, pvalue = pvalue, var_explained = var_explained))  # beta not included as using linear regressions where important value is R^2
  
}


# Define the lifestyle_info list
lifestyle_info <- list(
  list(field = 'p1558_i0', name = 'Alcohol intake frequency', type = 'categorical'),
  list(field = 'p20116_i0', name = 'Smoking status', type = 'categorical'),
  list(field = 'p23099_i0', name = '% body fat', type = 'continuous'),
  list(field = 'p1160_i0', name = 'Hours slept per day', type = 'continuous'),
  list(field = 'p22189', name = 'TDI', type = 'continuous'),
  list(field = 'p6138_i0', name = 'Education', type = 'categorical'),
  list(field = 'p884_i0', name = 'Days/week moderate activity', type = 'categorical'),
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

protein_vars_explained <- data.frame()

for (x in lifestyle_info){
  print(x$name)  # flush these to output log output
  flush.console()
  Sys.sleep(1)
  start_time <- Sys.time()
  model_result <- var_explained_by_covar(covars_df = covars,
                                         response_var_df = protein_data,
                                         covar = x$field,
                                         response_var = protein,
                                         covar_type = x$type,
                                         response_var_type = 'continuous')

  model_result$covar <- x$name
  
  model_result_df <- data.frame(model_result)
  
  print(Sys.time() - start_time)
  flush.console()
  Sys.sleep(1)
  
  # write_csv(model_result_df, paste0('lifestyle_effects/', disease, '_', x$field, '_lifestyle_effects_on_proteins.csv'))
  
  protein_vars_explained <- rbind(protein_vars_explained, model_result_df)
}

write_csv(protein_vars_explained, paste0('~/ch3_proteomics/3.1 covariates/lifestyle_effects/on_proteins/', protein, '_lifestyle_effects_on_proteins.csv'))


