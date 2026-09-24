# Minimally adjusted Firth regression for proteomics
library(logistf)
library(tidyverse)
library(lubridate)
library(progress)

# disease <- 'p131282'
disease_idx <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
# protein <- 'sf3b4'

# load in diseases
disease_info <- read.csv('~/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
diseases <- (disease_info %>%
               filter(disease_field != 'p131036'))$disease_field
disease <- diseases[disease_idx]
print(paste0(Sys.time(), ': Analysing proteins for ', disease))
flush.console()
Sys.sleep(1)

# load in unbalanced data
data <- read_csv('~/proteomics/proteomics_pan_ukbb_eur_post_regressions_robust_scaled_data.csv', num_threads = 1, show_col_types = FALSE)

# load in covariates
participant_covariates <- read_csv('~/proteomics/proteomics_pan_ukbb_eur_participant_covars.csv', num_threads = 1, show_col_types = FALSE)
protein_covariates <- read_csv('~/proteomics/proteomics_pan_ukbb_eur_protein_covars.csv', num_threads = 1, show_col_types = FALSE)

# load in phenotypes
phenotypes <- read_csv('~/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv', num_threads = 1, show_col_types = FALSE)

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
  inner_join(protein_covars_mean_dates, by = 'p30901_i0') %>%
  mutate(p53_i0_datetime = ymd(p53_i0),
         storage_time = as.numeric(difftime(mean_processing_start_date, p53_i0_datetime, units = "days")))

print(paste0(Sys.time(), ': Disease-inspecific data loaded.'))
flush.console()
Sys.sleep(1)

# disease-specific content
# load PCA data
pca <- read_table(paste0("/data/WHRI-Phenogenomics/krishna/pca/", disease, "/", disease, "_pca.eigenvec"), show_col_types = FALSE)

# Drop the '#FID' column and PC11 to PC20, then rename 'IID' to 'eid'
pca <- pca %>%
  select(-c(`#FID`, paste0("PC", 11:20))) %>%
  rename(eid = IID)

# add PCA data to data matrix
data_matrix <- data_matrix %>%
  inner_join(pca, by = 'eid')

print(paste0(Sys.time(), ': All data loaded. Fitting Firth models.'))
flush.console()
Sys.sleep(1)

all_firth_results <- data.frame()
proteins <- colnames(data)[2:ncol(data)]
bar <- progress_bar$new(
  format = "[:bar] :percent eta: :eta", 
  total = length(proteins),
  clear = FALSE, force = TRUE
)
for (protein in proteins){
  # remove samples with NA protein values
  data_matrix <- data_matrix %>%
    drop_na({{protein}})
  
  # set up model formula
  model_formula <- formula(paste0(disease, ' ~ ', protein, 
                                  ' + p21003_i0 + factor(p31) + factor(p54_i0) + storage_time + ',
                                  paste(paste0('PC', as.character(c(1:10))), collapse = ' + ')))
  
  # fit Firth model
  # max iterations for the model = 1000 (up from default of 25)
  # confidence intervals only calculated for the 2nd variable (the protein: Intercept is the 1st variable)
  firth_model <- logistf(model_formula, data_matrix, control = logistf.control(maxit = 1000), plconf = c(2))
  
  # export Firth model
  # saveRDS(firth_model, paste0('~/proteomics/association_analysis/firth_minimal/robjs/', protein, '_', disease, '.rds'))
  
  firth_results <- data.frame('disease' = disease, 'protein' = firth_model$terms, 'coefficient' = firth_model$coefficients, 'pval' = firth_model$prob) %>%
    filter(protein == {{protein}}) %>%
    mutate(beta = coefficient * sd(data_matrix %>% pull({{protein}})))
  
  all_firth_results <- rbind(all_firth_results, firth_results)
  
  bar$tick()
}

write_csv(all_firth_results, paste0('~/proteomics/association_analysis/firth_minimal/results/', disease, '_all_proteins.csv'))

print(paste0(Sys.time(), ': Analyses complete. Results and fitted models exported.'))
flush.console()
Sys.sleep(1)
