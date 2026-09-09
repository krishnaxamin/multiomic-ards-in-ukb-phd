# STEP 4: generate input info file
library(tidyverse)

# read in high-missingness samples from GWAS QC pipeline
high_missingness_samples <- read_tsv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/genomics/step4_high_missingness_samples_panukbb_eur.tsv',
                                     col_names = c('eid1', 'eid2', 'x'), col_select = any_of(c('eid1'))) %>%
  rename(eid = eid1)
# remove high-missingness samples from per-individual phenotype file
per_indiv_phenotypes <- read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv') %>%
  filter(!(eid %in% high_missingness_samples$eid))

# number of samples with GWAS actually done them
n_samples <- nrow(per_indiv_phenotypes)

# determine number of cases and controls for each disease, and set up info file
info_df <- data.frame('phenotype' = colnames(per_indiv_phenotypes %>% select(-eid)),
                      'cases' = colSums(per_indiv_phenotypes %>% select(-eid))) %>%
  mutate('controls' = n_samples - cases) %>%
  mutate(filename = paste0('~/ch2_genomics/2.11 local_genetic_correlations/int_data/gwas_summary_stats_for_lava/', phenotype, '.txt'))

write_tsv(info_df, '~/ch2_genomics/2.11 local_genetic_correlations/int_data/input.info.txt')
