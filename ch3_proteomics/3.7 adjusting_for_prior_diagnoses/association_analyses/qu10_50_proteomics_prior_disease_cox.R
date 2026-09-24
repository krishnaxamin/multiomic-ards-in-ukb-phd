# Lifestyle and prior disease-adjusted (with collinear variables removed) Cox regression for proteomics
library(Coxmos)
library(tidyverse)
library(lubridate)
library(progress)

# identify disease
# disease <- 'p131282'
disease <- as.character(commandArgs(trailingOnly = TRUE)[3])

# identify protein
# protein <- 'sf3b4'
protein_idx <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
proteins <- read.csv('~/data/internal/proteomics/protein_fields.csv')$protein_field
protein <- proteins[protein_idx]

# print to check disease-protein combo
print(paste0(Sys.time(), ': Analysing ', protein, ' for ', disease))
flush.console()
Sys.sleep(1)

protein_result_dir <- paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/cox_results/', disease)
protein_result_file <- paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/cox_results/', disease, '/', protein, '.csv')

# check if disease-protein results dir exists
if (!dir.exists(protein_result_dir)){
  dir.create(protein_result_dir)
} else {  # if it does exist, check if disease-protein already complete - if so, abort analysis
  if (file.exists(protein_result_file)) stop(paste0(disease, '-', protein, ' analysis already complete.'))
}

# load in unbalanced data
data <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_post_regressions_lifestyles_for_cox_robust_scaled_data.csv', num_threads = 1, show_col_types = FALSE)

# load in covariates
participant_covariates <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv', num_threads = 1, show_col_types = FALSE) %>%
  select(-p23099_i0)
protein_covariates <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv', num_threads = 1, show_col_types = FALSE)
lifestyles <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv', num_threads = 1, show_col_types = FALSE)

# load in phenotypes
phenotypes <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_future_phenotypes.csv', num_threads = 1, show_col_types = FALSE) %>%
  filter(disease == {{disease}})

# Convert the Processing_StartDate to a datetime
protein_covariates <- protein_covariates %>%
  mutate(processing_start_datetime = ymd(Processing_StartDate))

# Calculate the mean processing start date for each PlateID, then rename PlateID to 'p30901_i0'
protein_covars_mean_dates <- protein_covariates %>%
  group_by(PlateID) %>%
  summarise(mean_processing_start_date = mean(processing_start_datetime, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(p30901_i0 = PlateID)

# join all disease-inspecific data and convert the 'p53_i0' column to a datetime
data_matrix <- data %>%
  inner_join(participant_covariates, by = 'eid') %>%
  inner_join(phenotypes, by = 'eid') %>%
  inner_join(lifestyles, by = 'eid') %>%
  inner_join(protein_covars_mean_dates, by = 'p30901_i0') %>%
  mutate(p53_i0_datetime = ymd(p53_i0),
         storage_time = as.numeric(difftime(mean_processing_start_date, p53_i0_datetime, units = "days")))

print(paste0(Sys.time(), ': Disease-inspecific data loaded.'))
flush.console()
Sys.sleep(1)

# disease-specific content
# load genetic PCA data
pca <- read_table(paste0("~/data/internal/genomics/pca/", disease, "/", disease, "_pca.eigenvec"), show_col_types = FALSE)

# Drop the '#FID' column and PC11 to PC20, then rename 'IID' to 'eid'
pca <- pca %>%
  select(-c(`#FID`, paste0("PC", 11:20))) %>%
  rename(eid = IID)

# add PCA data to data matrix
data_matrix <- data_matrix %>%
  inner_join(pca, by = 'eid')

# load prior-disease PCA data
prior_disease_pcs <- read_csv(paste0('~/data/internal/proteomics/prior_disease_info/principal_components_before_ards/', disease, '.csv'), num_threads = 1, show_col_types = FALSE)
colnames(prior_disease_pcs) <- c('eid', paste0('disease_PC', as.character(c(1:(ncol(prior_disease_pcs) - 1)))))

# add prior-disease PCA data to data matrix
data_matrix <- data_matrix %>%
  inner_join(prior_disease_pcs, by = 'eid')

# specify which covariates are to be used, and select only those from the data_matrix
# note that this is the covariate set after collinearity analysis removed collinear variables
# also specify categorical covariates
covariates_for_cox <- c('p21003_i0', 'p31', 'p54_i0', 'p53_i0',
                        'p22189', 'p1438_i0', 'p1458_i0', 'p23099_i0',
                        paste0('p', as.character(seq(1289, 1319, 10)), '_i0'),
                        paste0('PC', as.character(c(1:10))),
                        paste0('disease_PC', as.character(c(1:33))))
columns_for_cox <- c('eid', 'time', 'event', covariates_for_cox, proteins)
categorical_covariates <- c('p31', 'p53_i0', 'p54_i0')

print(paste0(Sys.time(), ': All data loaded. Fitting Cox models.'))
flush.console()
Sys.sleep(1)

data_matrix_protein <- as.data.frame(data_matrix %>%
                                     drop_na({{protein}}) %>%
                                     select(all_of(c('eid', 'time', 'event', protein, covariates_for_cox))))

# make EIDs row names
rownames(data_matrix_protein) <- data_matrix_protein$eid
data_matrix_protein <- data_matrix_protein %>% select(!eid)

# isolate X and Y
x_data <- data_matrix_protein %>% select(all_of(c(protein, covariates_for_cox)))
y_data <- data_matrix_protein %>% select(time, event)

# calculate EPV for Cox
epv <- getEPV(x_data, y_data)

# different analyses based on EPV
if (epv >= 10) {  # normal, low-D Cox

  print(paste0(Sys.time(), ': EPV >= 10. Implementing normal Cox for ', protein, '.'))
  flush.console()
  Sys.sleep(0.1)
} else if (epv < 10) {  # high-D Cox -> sPLS-DRCOX
  
  print(paste0(Sys.time(), ': EPV < 10. High-D methods should be implemented for ', protein, ', but they are too time-inefficient for this use case. Implementing normal Cox instead.'))
  flush.console()
  Sys.sleep(0.1)

}

# turn categorical variables into factors
x_data <- x_data %>%
  mutate(across(all_of(categorical_covariates), as.factor))

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
    model_results_df <- data.frame(coefficient = numeric(),
                                   hazard_ratio = numeric(),
                                   pval = numeric(),
                                   protein = character(),
                                   epv = numeric(),
                                   disease = character())
  }
  
  # if Coxmos fails, but survival::coxph works, then results are contained in coxph_model from survival::coxph
  
  model_summary <- summary(coxph_model)
  print(paste0(protein, ' for ', disease, ': Coxmos failed, survival::coxph succeeded.'))
} else {
# if Coxmos has succeeded
  model_summary <- summary(cox_model$survival_model$fit)
}

# access results data
coefficient <- model_summary$coefficients[protein, 'coef']
hazard_ratio <- model_summary$coefficients[protein, 'exp(coef)']
pvalue <- model_summary$coefficients[protein, 'Pr(>|z|)']

model_result_df <- data.frame(coefficient = coefficient,
                              hazard_ratio = hazard_ratio,
                              pval = pvalue,
                              protein = protein,
                              epv = epv,
                              disease = disease) %>%
  mutate(beta = coefficient * sd(data_matrix_protein %>% pull({{protein}})))


# export per-disease Firth results
write_csv(model_result_df, protein_result_file)
