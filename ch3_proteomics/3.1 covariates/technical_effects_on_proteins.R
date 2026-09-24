' Getting covariate effects on proteomics data, using linear regression.
Script enabled for use on pre- and post-regressions data.'

library(tidyverse)

task_number <- as.integer(commandArgs(trailingOnly = TRUE)[1])

# load in specific protein data
proteins <- read.csv('~/data/internal/proteomics/protein_fields.csv')$protein_field
protein <- proteins[task_number]
protein_data <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_data.csv')[c('eid', protein)]

# load in participant covars
participant_covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv') %>%
  select(eid, p21003_i0, p31, p53_i0, p54_i0, p30901_i0, p30902_i0, Batch, plate_row, plate_column)

# load in protein covars
protein_covars <- read.csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv')

# adjust 'Assay' field in protein covars as required and isolate protein covars for the specific protein
if (protein %in% c('ervv_1', 'hla_a', 'hla_dra', 'hla_e')){
  single_protein_covars <- protein_covars %>%
    filter(tolower(Assay) == str_replace(protein, '_', '-'))
} else {
  single_protein_covars <- protein_covars %>%
    filter(tolower(Assay) == protein)
}

# use the PlateID to isolate the Processing_StartDate for the participants for this protein
eid_plate_protein_covars <- merge(participant_covars %>% select(eid, p30901_i0),
                                  single_protein_covars %>% select(-Assay, -UniProt),
                                  by.x = 'p30901_i0', by.y = 'PlateID')

# merge the Processing_StartDate info into the participant covars
proteomics_technical_covars <- merge(participant_covars,
                                     eid_plate_protein_covars %>% select(eid, Processing_StartDate),
                                     by = 'eid')

# calculate storage time
proteomics_technical_covars <- proteomics_technical_covars %>%
  mutate(storage_time = as.Date(Processing_StartDate) - as.Date(p53_i0))


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
  list(field = 'Processing_StartDate', name = 'date_processing', type = 'categorical'),
  list(field = 'storage_time', name = 'storage_time', type = 'continuous'),
  list(field = 'Batch', name = 'batch', type = 'categorical')
)
names(covariate_info) <- sapply(covariate_info, function(x) x$field)

protein_vars_explained <- data.frame()

for (x in covariate_info){
  print(x$name)  # flush these to output log output
  flush.console()
  Sys.sleep(1)
  start_time <- Sys.time()
  model_result <- var_explained_by_covar(covars_df = proteomics_technical_covars,
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
  
  # write_csv(model_result_df, paste0('covariate_effects/', disease, '_', x$field, '_covariate_effects_on_proteins.csv'))
  
  protein_vars_explained <- rbind(protein_vars_explained, model_result_df)
}

write_csv(protein_vars_explained, paste0('~/ch3_proteomics/3.1 covariates/technical_effects/on_proteins/', protein, '_technical_effects_on_proteins.csv'))
