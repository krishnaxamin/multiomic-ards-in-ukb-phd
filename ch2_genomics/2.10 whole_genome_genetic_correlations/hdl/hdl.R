# HDL

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("devtools")) install.packages("devtools")
library(devtools)
# install_github("zhenin/HDL/HDL")
library(HDL)

disease_pair <- commandArgs(trailingOnly = TRUE)[1]
disease1_field <- str_split(disease_pair, '-')[[1]][1]
disease2_field <- str_split(disease_pair, '-')[[1]][2]

print(disease1_field)
flush.console()
Sys.sleep(1)
print(disease2_field)
flush.console()
Sys.sleep(1)

# combine summary stats for all variants in disease1
disease1_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease1_field, '/', disease1_field, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1)
disease1_imputed_files <- list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease1_field, '/'),
                                     pattern = 'imputed.svat', full.names = TRUE)
disease1_imputed <- data.frame()
for (file in disease1_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease1_imputed <- rbind(disease1_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease1_data <- rbind(disease1_geno %>%
                         select(MarkerID, Allele1, Allele2, BETA, SE, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, b = BETA, se = SE, SNP = MarkerID),  # Allele2 is the allele for which stats are calculated, A1 = effect allele
                       disease1_imputed %>%
                         select(MarkerID, Allele1, Allele2, BETA, SE, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, b = BETA, se = SE, SNP = MarkerID)) %>%  # Allele2 is the allele for which stats are calculated, A1 = effect allele
  mutate(N = N_case + N_ctrl) %>%
  select(-N_case, -N_ctrl) %>%
  distinct(SNP, A1, A2, .keep_all = TRUE) %>%  # this is slightly redundant because HDL doesn't preserve multiple alleles at a locus anyway
  select(SNP, A1, A2, N, b, se)

print(summary(disease1_data))

# combine summary stats for all variants in disease2
disease2_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/', disease2_field, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1)
disease2_imputed_files <- list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/'),
                                     pattern = 'imputed.svat', full.names = TRUE)
disease2_imputed <- data.frame()
for (file in disease2_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease2_imputed <- rbind(disease2_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease2_data <- rbind(disease2_geno %>%
                         select(MarkerID, Allele1, Allele2, BETA, SE, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, b = BETA, se = SE, SNP = MarkerID),  # Allele2 is the allele for which stats are calculated, A1 = effect allele
                       disease2_imputed %>%
                         select(MarkerID, Allele1, Allele2, BETA, SE, N_case, N_ctrl) %>%
                         rename(A1 = Allele2, A2 = Allele1, b = BETA, se = SE, SNP = MarkerID)) %>%  # Allele2 is the allele for which stats are calculated, A1 = effect allele
  mutate(N = N_case + N_ctrl) %>%
  select(-N_case, -N_ctrl) %>%
  distinct(SNP, A1, A2, .keep_all = TRUE) %>%  # this is slightly redundant because HDL doesn't preserve multiple alleles at a locus anyway
  select(SNP, A1, A2, N, b, se)

print(summary(disease2_data))

hdl_test <- HDL.rg(disease1_data, disease2_data, 
                   LD.path = '~/data/external/genomics/hdl/UKB_imputed_SVD_eigen99_extraction',
                   intercept.output = TRUE,
                   output.file = paste0('~/ch2_genomics/2.10 whole_genome_genetic_correlations/hdl/out_files/', disease1_field, '-', disease2_field, '.hdl.out'))
saveRDS(hdl_test, paste0('~/ch2_genomics/2.10 whole_genome_genetic_correlations/hdl/results/', disease1_field, '-', disease2_field, '.hdl'))
