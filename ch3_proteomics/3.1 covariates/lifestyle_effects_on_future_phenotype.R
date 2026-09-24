' Getting lifestyle effects on proteomics future phenotypes, using Cox regressions '

library(tidyverse)
library(Coxmos)

covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv')


# Define the lifestyle_info list
lifestyle_info <- list(
  list(field = 'p1558_i0', name = 'Alcohol intake frequency', type = 'categorical'),
  list(field = 'p20116_i0', name = 'Smoking status', type = 'categorical'),
  list(field = 'p23099_i0', name = '% body fat', type = 'continuous'),
  list(field = 'p1160_i0', name = 'Hours slept per day', type = 'continuous'),
  list(field = 'p22189', name = 'TDI', type = 'continuous'),
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

# read in disease data
qu10_50_disease_fields <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')$disease_field
phenotype <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_future_phenotypes.csv')[c('eid', disease)]

phenotype_vars_explained <- data.frame()
for (disease in qu10_50_disease_fields) {
  print(disease)  # flush these to output log output
  flush.console()
  Sys.sleep(0.1)
  
  # select phenotype for this one disease
  disease_phenotype <- phenotype  %>%
    filter(disease == {{disease}}) %>%
    select(eid, time, event)
  
  # merge phenotype and covariates
  data_df <- merge(disease_phenotype, covars)
  
  # make EIDs row names
  rownames(data_df) <- data_df$eid
  data_df <- data_df %>% select(!eid)
  
  if (var(phenotype['time']) != 0) {
    for (x in lifestyle_info) {
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
}

write_csv(phenotype_vars_explained, '~/ch3_proteomics/3.1 covariates/lifestyle_effects/lifestyle_effects_on_future_phenotype.csv')
