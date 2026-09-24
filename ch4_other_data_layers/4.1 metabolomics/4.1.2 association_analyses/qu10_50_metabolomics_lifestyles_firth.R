# Lifestyle-adjusted (with collinear variables removed) Firth regression for metabolomics
library(logistf)
library(tidyverse)
library(lubridate)
library(progress)

# disease <- 'p131282'
disease_idx <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
# metabolite <- 'sf3b4'

# load in diseases
disease_info <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
diseases <- disease_info$disease_field
disease <- diseases[disease_idx]
print(paste0(Sys.time(), ': Analysing metabolites for ', disease))
flush.console()
Sys.sleep(1)

# load in unbalanced data
data <- read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_post_regressions_robust_scaled_data.csv', num_threads = 1, show_col_types = FALSE)

# load in covariates
covariates <- read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_covars.csv', num_threads = 1, show_col_types = FALSE) %>%
  select(-p23099_i0)
lifestyles <- read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_lifestyles.csv', num_threads = 1, show_col_types = FALSE)

# load in phenotypes
phenotypes <- read_csv('~/data/internal/metabolomics/metabolomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv', num_threads = 1, show_col_types = FALSE)

# join all disease-inspecific data and convert the 'p53_i0' column to a datetime
data_matrix <- data %>%
  inner_join(covariates, by = 'eid') %>%
  inner_join(phenotypes, by = 'eid') %>%
  inner_join(lifestyles, by = 'eid')

print(paste0(Sys.time(), ': Disease-inspecific data loaded.'))
flush.console()
Sys.sleep(1)

# disease-specific content
# load PCA data
pca <- read_table(paste0("~/data/internal/genomics/pca/", disease, "/", disease, "_pca.eigenvec"), show_col_types = FALSE)

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
metabolites <- colnames(data)[2:ncol(data)]
bar <- progress_bar$new(
  format = "[:bar] :percent eta: :eta", 
  total = length(metabolites),
  clear = FALSE, force = TRUE
)
for (metabolite in metabolites){
  # remove samples with NA metabolite values
  data_matrix <- data_matrix %>%
    drop_na({{metabolite}})
  
  # set up model formula - collinear factors removed (CHECK)
  model_formula <- formula(paste0(disease, ' ~ ', metabolite, 
                                  ' + p21003_i0 + factor(p31) + factor(p54_i0) + storage_time + p22189 + p884_i0 +
                                  p1438_i0 + p1458_i0 + p1528_i0 + p23099_i0 + ',
                                  paste(paste0('p', as.character(seq(1289, 1319, 10)), '_i0'), collapse = ' + '), ' + ',
                                  paste(paste0('PC', as.character(c(1:10))), collapse = ' + ')))
  
  # fit Firth model
  # max iterations for the model = 1000 (up from default of 25)
  # confidence intervals only calculated for the 2nd variable (the metabolite: Intercept is the 1st variable)
  firth_model <- logistf(model_formula, data_matrix, control = logistf.control(maxit = 1000), plconf = c(2))
  
  firth_results <- data.frame('disease' = disease, 'metabolite' = firth_model$terms, 'coefficient' = firth_model$coefficients, 'pval' = firth_model$prob) %>%
    filter(metabolite == {{metabolite}}) %>%
    mutate(beta = coefficient * sd(data_matrix %>% pull({{metabolite}})))
  
  all_firth_results <- rbind(all_firth_results, firth_results)
  
  bar$tick()
}

write_csv(all_firth_results, paste0('~/ch4_other_data_layers/4.1 metabolomics/4.1.2 association_analyses/firth_results/', disease, '_all_metabolites.csv'), num_threads = 1)

print(paste0(Sys.time(), ': Analyses complete. Results and fitted models exported.'))
flush.console()
Sys.sleep(1)
