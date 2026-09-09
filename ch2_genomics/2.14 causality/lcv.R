# LCV

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

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
disease1_imputed_files <- list.files(path = paste0('ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease1_field, '/'), 
                                     pattern = 'imputed.svat', full.names = TRUE)
disease1_imputed <- data.frame()
for (file in disease1_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease1_imputed <- rbind(disease1_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease1_data <- rbind(disease1_geno %>%
                         filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354)) %>%  # remove MHC
                         filter(abs(AF_Allele2 - 0.5) < 0.45) %>%  # MAF > 5%
                         select(CHR, POS, Allele1, Allele2, BETA, SE), # Allele2 is the allele for which stats are calculated, A1 = effect allele
                       disease1_imputed %>%
                         filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354)) %>%  # remove MHC
                         filter(abs(AF_Allele2 - 0.5) < 0.45) %>%  # MAF > 5%
                         select(CHR, POS, Allele1, Allele2, BETA, SE)) %>%  # Allele2 is the allele for which stats are calculated, A1 = effect allele
  mutate(z_score = BETA / SE) %>%
  select(-BETA, -SE) %>%
  rowwise() %>%
  mutate(SNP = paste(CHR, POS, paste(sort(c(Allele1, Allele2)), collapse = ':'), sep = ':')) %>% 
  ungroup() %>%
  distinct(SNP, .keep_all = TRUE) %>%  # this is slightly redundant because HDL doesn't preserve multiple alleles at a locus anyway
  select(SNP, CHR, POS, z_score) %>%
  rename(z_score_disease1 = z_score)

print(summary(disease1_data))

# combine summary stats for all variants in disease2
disease2_geno <- read_tsv(paste0('ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/', disease2_field, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1)
disease2_imputed_files <- list.files(path = paste0('ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/'), 
                                     pattern = 'imputed.svat', full.names = TRUE)
disease2_imputed <- data.frame()
for (file in disease2_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease2_imputed <- rbind(disease2_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease2_data <- rbind(disease2_geno %>%
                         filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354)) %>%  # remove MHC
                         filter(abs(AF_Allele2 - 0.5) < 0.45) %>%  # MAF > 5%
                         select(CHR, POS, Allele1, Allele2, BETA, SE),  # Allele2 is the allele for which stats are calculated, A1 = effect allele
                       disease2_imputed %>%
                         filter(!(CHR == 6 & POS >= 28477797 & POS <= 33448354)) %>%  # remove MHC
                         filter(abs(AF_Allele2 - 0.5) < 0.45) %>%  # MAF > 5%
                         select(CHR, POS, Allele1, Allele2, BETA, SE)) %>%  # Allele2 is the allele for which stats are calculated, A1 = effect allele
  mutate(z_score = BETA / SE) %>%
  select(-BETA, -SE) %>%
  rowwise() %>%
  mutate(SNP = paste(CHR, POS, paste(sort(c(Allele1, Allele2)), collapse = ':'), sep = ':')) %>% 
  ungroup() %>%
  distinct(SNP, .keep_all = TRUE) %>%  # this is slightly redundant because HDL doesn't preserve multiple alleles at a locus anyway
  select(SNP, CHR, POS, z_score) %>%
  rename(z_score_disease2 = z_score)

print(summary(disease2_data))

# read in LD data
ld <- read_tsv('~/data/external/UKBB.EUR.l2.ldscore', show_col_types = FALSE, num_threads = 1)
ld_cleaned <- ld %>%
  separate(SNP, into = c("CHR", "POS", "a1", "a2"), sep = ":") %>%
  rowwise() %>%
  mutate(
    alleles = list(sort(c(a1, a2))),
    a1 = alleles[1],
    a2 = alleles[2],
    SNP = paste(CHR, POS, a1, a2, sep = ":")
  ) %>%
  select(SNP, CHR, POS, L2) %>%
  ungroup()

# merge disease1 data, disease2 data, and LD data
data <- disease1_data %>%
  mutate(across(c(CHR, POS), as.numeric)) %>%
  inner_join(disease2_data %>% mutate(across(c(CHR, POS), as.numeric))) %>%
  inner_join(ld_cleaned %>% mutate(across(c(CHR, POS), as.numeric))) %>%
  arrange(CHR, POS)

# set up LCV
setwd("~/ch2_genomics/2.14 causality/oconnor")
source("~/ch2_genomics/2.14 causality/oconnor/RunLCV.R")

# run LCV and export
lcv = RunLCV(data$L2, data$z_score_disease1, data$z_score_disease2)
saveRDS(lcv, file = paste0('~/ch2_genomics/2.14 causality/results/', disease1_field, '_',
                           disease2_field, '_lcv.rds'))
