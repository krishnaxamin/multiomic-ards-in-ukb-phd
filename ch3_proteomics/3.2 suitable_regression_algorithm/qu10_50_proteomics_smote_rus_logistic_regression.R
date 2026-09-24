' Performing logistic regression on SMOTE-RUS balanced and robust-scaled data. 
Implements penalised elastic-net regression if alpha < 0.05 or model warnings arise.'

library(tidyverse)
library(glmnet)

# disease_protein <- commandArgs(trailingOnly = TRUE)[1]
# disease <- str_split(disease_protein, '-')[[1]][1]
protein <- commandArgs(trailingOnly = TRUE)[1]

qu10_50_disease_fields <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')$disease_field

available_balanced_data <- list.files(path = '~/ch3_proteomics/3.2 suitable_regression_algorithm/balanced_data/', pattern = paste0(protein, '_balanced_data.csv'))

results_df <- data.frame()
for (disease in qu10_50_disease_fields){
  if (paste0(disease, '-', protein, '_balanced_data.csv') %in% available_balanced_data){
    balanced_data <- read.csv(paste0('~/ch3_proteomics/3.2 suitable_regression_algorithm/balanced_data/', disease, '-', protein, '_balanced_data.csv'))
    
    model_formula <- formula(paste0(disease, ' ~ ', protein, ' + p21003_i0 + factor(p31) + factor(p54_i0) + storage_time + ',paste(paste0('PC', as.character(c(1:10))), collapse = ' + ')))
    
    start_time <- Sys.time()
    
    # catch Warning messages if they arise
    warning_counter = 0
    tryCatch(expr = {model <- glm(formula = model_formula, data = balanced_data, family = 'binomial')}, 
             warning = function(w){
               warning_counter = warning_counter + 1
             }) 
    suppressWarnings(model <- glm(formula = model_formula, data = balanced_data, family = 'binomial'))
    
    # get stats from normal logistic regression
    model_summary <- summary(model)
    coefficient <- model_summary$coefficients[, 'Estimate'][2]
    beta <- coefficient * sd(balanced_data[, protein])  # beta is response-standardised only, as (1) s.d. of disease state always the same (2) what does the s.d. of a binary variable mean when interpreting beta
    pvalue <- model_summary$coefficients[, 'Pr(>|z|)'][2]
    
    # penalised regression if required
    penalised_model_coeff <- NA
    penalised_beta <- NA
    best_lambda <- NA
    if ((warning_counter == 1) | (pvalue < (0.05/ (67 * 2920)))){
      penalised_data <- balanced_data
      
      # turn categorical variables into factors and set the reference factor
      penalised_data$p31 <- relevel(as.factor(penalised_data$p31), ref = 'Female')
      penalised_data$p54_i0 <- relevel(as.factor(penalised_data$p54_i0), ref = 'Barts')
      
      # set up model matrix
      penalised_model_formula <- formula(paste0(' ~ ', protein, ' + p21003_i0 + factor(p31) + factor(p54_i0) + storage_time + ', paste(paste0('PC', as.character(c(1:10))), collapse = ' + '), ' - 1'))
      model_matrix <- model.matrix(model_formula, penalised_data)
      
      # get out phenotype
      phenotype <- balanced_data %>%
        select({{disease}})
      
      # glmnet cv
      glmnet_cv <- cv.glmnet(x = model_matrix, y = as.matrix(phenotype), family = 'binomial', alpha = 0.5)
      
      # identify the best lambda and isolate coefficients from using that lambda
      best_lambda <- glmnet_cv$lambda.min
      glmnet_coefs <- coef(glmnet_cv, s = 'lambda.min')
      
      # obtain protein coefficent and beta
      if (protein %in% glmnet_coefs@Dimnames[[1]][glmnet_coefs@i]){
        # get the protein's coefficient from the penalised model result
        protein_coeff_idx <- match(match(protein, glmnet_coefs@Dimnames[[1]]), glmnet_coefs@i)
        penalised_model_coeff <- glmnet_coefs@x[protein_coeff_idx]
        # standardise to beta
        penalised_beta <- penalised_model_coeff * sd(balanced_data[, protein])  # beta is response-standardised only, as (1) s.d. of disease state always the same (2) what does the s.d. of a binary variable mean when interpreting beta
      }
    }
    
    # print time taken
    print(Sys.time() - start_time)
    flush.console()
    Sys.sleep(1)
    
    result_df <- data.frame(disease_field = disease,
                            protein_field = protein,
                            log_odds = coefficient,
                            beta = beta,
                            pvalue = pvalue,
                            penalised_log_odds = penalised_model_coeff,
                            penalised_beta = penalised_beta,
                            penalised_opt_lambda = best_lambda)
    
    results_df <- rbind(results_df, result_df)
  } else {
    print(paste0('Balanced data for ', disease, ' is not available, because there are 0-1 prevalent cases in this cohort.'))
  }
}
write_csv(results_df, paste0('~/ch3_proteomics/3.2 suitable_regression_algorithm/association_results/', protein, '_unpenalised_elastic_net_results.csv'))
