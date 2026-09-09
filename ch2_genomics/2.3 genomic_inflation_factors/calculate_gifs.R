# Calculate genomic inflation factor (GIF) for SAIGE and REGENIE trial runs.

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

chr_list <- c(1:22, 'XY')

disease <- commandArgs(trailingOnly = TRUE)[1]

print(disease)

saige_geno_complete <- data.frame()
saige_imputed_complete <- data.frame()
regenie_imputed_complete <- data.frame()

# collate per-chromosome results (variant ID, alleles, p_val, test stat) for SAIGE and REGENIE
for (chr in chr_list){
  saige_geno <- read.csv(paste0('~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/saige_svat/',
                                     disease, '/', disease, '_c', chr, '_pan-ukbb-eur.geno.svat'), sep = '') %>%
    rowwise() %>%
    mutate(allele_pair = paste(sort(c(Allele1, Allele2)), collapse = '/')) %>%
    mutate(algo = 'SAIGE') %>%
    select(MarkerID, allele_pair, p.value, algo, Tstat)
  saige_geno_complete <- rbind(saige_geno_complete, saige_geno)
  
  saige_imputed <- read.csv(paste0('~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/saige_svat/',
                                        disease, '/', disease, '_c', chr, '_pan-ukbb-eur.imputed.svat'), sep = '') %>%
    filter(abs(AF_Allele2) - 0.5 < 0.49) %>%
    rowwise() %>%
    mutate(allele_pair = paste(sort(c(Allele1, Allele2)), collapse = '/')) %>%
    mutate(algo = 'SAIGE') %>%
    select(MarkerID, allele_pair, p.value, algo, Tstat)
  saige_imputed_complete <- rbind(saige_imputed_complete, saige_imputed)
  
  regenie_imputed <- read.csv(paste0('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step2/',
                                     disease, '/', disease, '_c', chr, '_pan-ukbb-eur.imputed.step2_', disease, '.regenie'), sep = '') %>%
    filter(abs(A1FREQ) - 0.5 < 0.49) %>%
    rowwise() %>%
    mutate(allele_pair = paste(sort(c(ALLELE0, ALLELE1)), collapse = '/')) %>%
    mutate(algo = 'REGENIE') %>%
    select(ID, allele_pair, LOG10P, algo, CHISQ)
  regenie_imputed_complete <- rbind(regenie_imputed_complete, regenie_imputed)
}

# genotyped results from REGENIE run on all chromosomes at once
regenie_geno_complete <- read.csv(paste0('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/regenie_step2/',
                                         disease, '/', disease, '_pan-ukbb-eur.geno.step2_', disease, '.regenie'), sep = '') %>%
  rowwise() %>%
  mutate(allele_pair = paste(sort(c(ALLELE0, ALLELE1)), collapse = '/')) %>%
  mutate(algo = 'REGENIE') %>%
  select(ID, allele_pair, LOG10P, algo, CHISQ)

# collate SAIGE genotyped and imputed results and process to get obs and exp log10P for QQ plot
all_saige_vars <- rbind(saige_geno_complete, saige_imputed_complete) %>% 
  distinct(MarkerID, allele_pair, .keep_all = TRUE) 
all_saige_vars <- all_saige_vars %>%
  select(p.value, algo, Tstat) %>%
  arrange(p.value) %>%
  mutate(obs_log = -1 * log10(p.value))
all_saige_vars$expected = 1:nrow(all_saige_vars)
all_saige_vars <- all_saige_vars %>%
  mutate(expected_log = -1 * log10(expected / (nrow(all_saige_vars) + 1))) %>%
  mutate(chisq_from_pval = qchisq(1 - p.value, 1))

# write_csv(all_saige_vars, 'gwas_comparisons/all_saige_vars_test2.csv')

# collate REGENIE genotyped and imputed results and process to get obs and exp log10P for QQ plot
all_regenie_vars <- rbind(regenie_geno_complete, regenie_imputed_complete) %>% 
  distinct(ID, allele_pair, .keep_all = TRUE) 
all_regenie_vars <- all_regenie_vars %>%
  select(LOG10P, algo, CHISQ) %>%
  arrange(desc(LOG10P))
all_regenie_vars$expected = 1:nrow(all_regenie_vars)
all_regenie_vars <- all_regenie_vars %>%
  mutate(expected_log = -1 * log10(expected / (nrow(all_regenie_vars) + 1))) %>%
  rename(obs_log = LOG10P) %>%
  mutate(chisq_from_pval = qchisq(1 - 10 ** (-1 * obs_log), 1))

# calculate Genomic Inflation Factors and export

saige_gif_pval <- median(all_saige_vars$chisq_from_pval) / qchisq(0.5, 1)

regenie_gif_chisq <- median(all_regenie_vars$CHISQ) / qchisq(0.5, 1)
regenie_gif_pval <- median(all_regenie_vars$chisq_from_pval) / qchisq(0.5, 1)

write_csv(data.frame(saige_gif_pval = saige_gif_pval,
                     regenie_gif_chisq = regenie_gif_chisq,
                     regenie_gif_pval = regenie_gif_pval), 
          paste0('~/ch2_genomics/2.3 genomic_inflation_factors/disease_gifs/', disease, '_gifs.csv'))
