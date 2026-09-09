# STEP 1: sumstats -> LAVA sumstats, per disease
# SNP; A1 (effect); A2 (reference); N (n_samples); B (effect size); P (pval)

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# digested from ~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt
disease <- commandArgs(trailingOnly = TRUE)[1]

print(disease)
flush.console()
Sys.sleep(1)

# combine summary stats for all variants in disease
cols_to_read <- c('MarkerID', 'Allele1', 'Allele2', 'BETA', 'p.value', 'N_case', 'N_ctrl')
disease_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease, '/', disease, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1, col_select = any_of(cols_to_read))
disease_imputed_files <- list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease, '/'),
                                     pattern = 'imputed.svat', full.names = TRUE)
disease_imputed <- data.frame()
for (file in disease_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease_imputed <- rbind(disease_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1, 
                                                       col_select = any_of(c('imputationInfo', cols_to_read))) %>% 
                               filter(imputationInfo > 0.8) %>%
                               select(-imputationInfo))
  }
}
disease_data <- rbind(disease_geno %>%
                         select(MarkerID, Allele1, Allele2, BETA, p.value, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, B = BETA, P = p.value, SNP = MarkerID),  # Allele2 is the allele for which stats are calculated, A1 = effect allele
                       disease_imputed %>%
                         select(MarkerID, Allele1, Allele2, BETA, p.value, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, B = BETA, P = p.value, SNP = MarkerID)) %>%  # Allele2 is the allele for which stats are calculated, A1 = effect allele
  mutate(N = N_case + N_ctrl) %>%
  select(-N_case, -N_ctrl) %>%
  distinct(SNP, A1, A2, .keep_all = TRUE) %>% 
  select(SNP, A1, A2, N, B, P)

write_tsv(disease_data, paste0('~/ch2_genomics/2.11 local_genetic_correlations/int_data/gwas_summary_stats_for_lava/', disease, '.txt'))
