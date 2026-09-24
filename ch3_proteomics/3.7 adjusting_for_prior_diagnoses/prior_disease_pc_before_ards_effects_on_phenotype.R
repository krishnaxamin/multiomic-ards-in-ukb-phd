' Getting lifestyle effects on proteomics future phenotypes, using linear regression '

library(tidyverse)
library(Coxmos)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

qu10_50_disease_fields <- read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt', num_threads = 1, show_col_types = FALSE)$disease_field
disease <- qu10_50_disease_fields[task_number]
phenotype <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_future_phenotypes.csv', num_threads = 1, show_col_types = FALSE)
pcs <- read_csv(paste0('~/data/internal/proteomics/prior_disease_info/principal_components_before_ards/', disease, '.csv'), num_threads = 1, show_col_types = FALSE)

# Define the pcs_info list
pcs_info <- lapply(1:(ncol(pcs) - 1), function(x) list(field = paste0('PC', x), name = paste0('PC', x), type = 'continuous'))
names(pcs_info) <- sapply(pcs_info, function(x) x$field)

print(disease)  # flush these to output log output
flush.console()
Sys.sleep(0.1)

# select phenotype for this one disease
disease_phenotype <- phenotype  %>%
  filter(disease == {{disease}}) %>%
  select(eid, time, event)

# merge phenotype and covariates
data_df <- merge(disease_phenotype, pcs)

# make EIDs row names
rownames(data_df) <- data_df$eid
data_df <- data_df %>% select(!eid)

phenotype_vars_explained <- data.frame()
if (var(phenotype['time']) != 0) {
  for (x in pcs_info) {
    # if (x$field != 'storage_time'){
    #   next
    # }
    print(x$name)  # flush these to output log output
    flush.console()
    Sys.sleep(0.1)
    start_time <- Sys.time()
    
    # isolate X and Y
    x_data <- data_df[, x$field, drop = FALSE]
    y_data <- data_df %>% select(time, event)
    
    # turn categorical variables into factors
    if (x$type == 'categorical'){
      x_data[, x$field] <- as.factor(x_data[, x$field])
    }
    
    # run Cox
    cox_model <- tryCatch(
      expr = {coxmos(method = 'cox', 
                     X = x_data,
                     Y = y_data)},
      error = function(e){
        message(paste0("COX: ", e))
        return(NULL)
      }
    )
    # cox_model <- coxmos(method = 'cox', 
    #                     X = x_data,
    #                     Y = y_data)
    if (is.null(cox_model)){
      
      y_time_scaled <- scale(y_data[,"time", drop = FALSE], scale = FALSE)
      y_data_scaled <- y_data
      y_data_scaled[, 'time'] <- y_time_scaled
      
      if (x$type == 'continuous'){
        x_data_scaled <- scale(x_data, scale = FALSE)
      } else {
        x_data_scaled <- x
      }
      
      survival_d <- as.data.frame(cbind(x_data_scaled, y_data_scaled))
      
      coxph_model <- tryCatch(
        expr = {
          survival::coxph(formula = survival::Surv(time,event) ~ .,
                          data = survival_d,
                          ties = "efron",
                          singular.ok = TRUE,
                          robust = TRUE,
                          nocenter = rep(1, ncol(X)),
                          model = TRUE, x = TRUE)
        },
        error = function(e){
          message(paste0("COX: ", e))
          return(NULL)
        } 
      )
    }
    
    if (is.null(cox_model)){
      # if Coxmos fails
      
      if (is.null(coxph_model)){
        # if survival::coxph fails, then print failure and skip to next 
        print(paste0(protein, ' for ', disease, ': both Cox failed.'))
        flush.console()
        Sys.sleep(0.1)
        next
      }
      
      # if Coxmos fails, but survival::coxph works, then results are contained in coxph_model from survival::coxph
      
      model_summary <- summary(coxph_model)
      print(paste0(protein, ' for ', disease, ': Coxmos failed, survival::coxph succeeded.'))
    } else {
      # if Coxmos has succeeded
      model_summary <- summary(cox_model$survival_model$fit)
    }
    
    
    coefficient <- model_summary$coefficients[x$field, 'coef']
    hazard_ratio <- model_summary$coefficients[x$field, 'exp(coef)'] 
    pvalue <- model_summary$coefficients[x$field, 'Pr(>|z|)']
    
    model_result_df <- data.frame(coefficient = coefficient,
                                  hazard_ratio = hazard_ratio,
                                  pvalue = pvalue,
                                  covar = x$name,
                                  disease = disease)
    
    print(Sys.time() - start_time)
    flush.console()
    Sys.sleep(0.1)
    
    phenotype_vars_explained <- rbind(phenotype_vars_explained, model_result_df)
  }
}

write_csv(phenotype_vars_explained, paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/before_ards_effect_on_phenotype/', disease, '_prior_disease_pc_effects_on_phenotype.csv'))
