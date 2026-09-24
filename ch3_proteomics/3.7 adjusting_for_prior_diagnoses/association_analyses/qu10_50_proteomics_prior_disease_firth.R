# Lifestyle and prior disease-adjusted (with collinear variables removed) Firth regression for proteomics
library(logistf)
library(tidyverse)
library(lubridate)
library(progress)

# identify disease
# disease <- 'p131282'
disease <- as.character(commandArgs(trailingOnly = TRUE)[3])
if (disease == 'p131036') stop('p131036 has no prevalent cases.')

# identify protein
# protein <- 'sf3b4'
protein_idx <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
proteins <- read.csv('~/data/internal/proteomics/protein_fields.csv')$protein_field
protein <- proteins[protein_idx]

# print to check disease-protein combo
print(paste0(Sys.time(), ': Analysing ', protein, ' for ', disease))
flush.console()
Sys.sleep(1)

protein_result_dir <- paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/firth_results/', disease)
protein_result_file <- paste0('~/ch3_proteomics/3.7 adjusting_for_prior_diagnoses/association_analyses/firth_results/', disease, '/', protein, '.csv')

# check if disease-protein results dir exists
if (!dir.exists(protein_result_dir)){
  dir.create(protein_result_dir)
} else {  # if it does exist, check if disease-protein already complete - if so, abort analysis
  if (file.exists(protein_result_file)) stop(paste0(disease, '-', protein, ' analysis already complete.'))
}

# load in unbalanced data
data <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_post_regressions_lifestyles_robust_scaled_data.csv', num_threads = 1, show_col_types = FALSE)

# load in covariates
participant_covariates <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv', num_threads = 1, show_col_types = FALSE) %>%
  select(-p23099_i0)
protein_covariates <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv', num_threads = 1, show_col_types = FALSE)
lifestyles <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_lifestyles.csv', num_threads = 1, show_col_types = FALSE)

# load in phenotypes
phenotypes <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv', num_threads = 1, show_col_types = FALSE)

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
prior_disease_pcs <- read_csv(paste0('~/data/internal/proteomics/prior_disease_info/principal_components_before_blood_before_ards/', disease, '.csv'), num_threads = 1, show_col_types = FALSE)
colnames(prior_disease_pcs) <- c('eid', paste0('disease_PC', as.character(c(1:(ncol(prior_disease_pcs) - 1)))))

# add prior-disease PCA data to data matrix
data_matrix <- data_matrix %>%
  inner_join(prior_disease_pcs, by = 'eid')

print(paste0(Sys.time(), ': All data loaded. Fitting Firth models.'))
flush.console()
Sys.sleep(1)

# remove samples with NA protein values
data_matrix <- data_matrix %>%
  drop_na({{protein}})

# set up model formula - collinear factors removed (CHECK)
model_formula <- formula(paste0(disease, ' ~ ', protein,
                                ' + p21003_i0 + factor(p31) + factor(p54_i0) + storage_time +
                                p22189 + p884_i0 + p1438_i0 + p1458_i0 + p23099_i0 + ',
                                paste(paste0('p', as.character(seq(1309, 1319, 10)), '_i0'), collapse = ' + '), ' + ',
                                paste(paste0('PC', as.character(c(1:10))), collapse = ' + '), ' + ',
                                paste(paste0('disease_PC', as.character(c(1:(ncol(prior_disease_pcs) - 1)))), collapse = ' + ')))

# fit Firth model - skipping if a 'determinant is 0' error is thrown and recording which disease-protein combo threw the error
# max iterations for the model = 1000 (up from default of 25)
# confidence intervals only calculated for the 2nd variable (the protein: Intercept is the 1st variable)
# ensures firth_model is always storing something
firth_model <- tryCatch(expr = {logistf(model_formula, data_matrix, control = logistf.control(maxit = 1000), plconf = c(2))},
         error = function(e){
           print(paste0(protein, ' threw: ', e))
           flush.console()
           Sys.sleep(1)
           return(NULL)  # makes firth_model NULL if error thrown
         })
if (is.null(firth_model)){
  write_csv(data.frame(disease=character(), protein=character(), coefficient=numeric(), pval=numeric()), protein_result_file)
  stop()
}

firth_results <- data.frame('disease' = disease, 'protein' = firth_model$terms, 'coefficient' = firth_model$coefficients, 'pval' = firth_model$prob) %>%
  filter(protein == {{protein}}) %>%
  mutate(beta = coefficient * sd(data_matrix %>% pull({{protein}})))

# export per-disease Firth results
write_csv(firth_results, protein_result_file)
