' Getting covariate effects on proteomics future phenotypes, using Cox regression '

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
  mutate(storage_time = as.integer(as.Date(mean_processing_start_date) - as.Date(p53_i0)))


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

# read in disease data
qu10_50_disease_fields <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')$disease_field
phenotype <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_future_phenotypes.csv')

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
    for (x in covariate_info) {
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
      
      if (is.null(coxph_model)){
        # if survival::coxph fails, then print failure and skip to next 
        print(paste0(x$name, ' for ', disease, ': both Cox failed.'))
        flush.console()
        Sys.sleep(0.1)
        next
      }
      
      if (is.null(cox_model)){
        # if Coxmos fails but survival::coxph succeeds
        model_summary <- summary(coxph_model)
        print(paste0(x$name, ' for ', disease, ': Coxmos failed, survival::coxph succeeded.'))
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

write_csv(phenotype_vars_explained, '~/ch3_proteomics/3.1 covariates/technical_effects/technical_effects_on_future_phenotype.csv')
