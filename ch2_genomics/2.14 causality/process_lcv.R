# process LCV results

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

disease_pairs <- read.csv('./disease_pairs.txt')$disease_pair

# collate
lcv_results <- data.frame()
for (results_file in list.files(paste0('/ch2_genomics/2.14 causality/results/'))){
  disease1 <- str_split(results_file, '_')[[1]][1]
  disease2 <- str_split(results_file, '_')[[1]][2]
  
  lcv <- readRDS(paste0('~/ch2_genomics/2.14 causality/results/', results_file))
  
  lcv_results <- rbind(lcv_results,
                       data.frame('disease1' = disease1,
                                  'disease2' = disease2,
                                  'zscore' = lcv$zscore,
                                  'gcp0_pval' = lcv$pval.gcpzero.2tailed,
                                  'gcp_est' = lcv$gcp.pm,
                                  'gcp_se' = lcv$gcp.pse,
                                  'corr_est' = lcv$rho.est,
                                  'corr_err' = lcv$rho.err,
                                  'gcppos1_pval' = lcv$pval.fullycausal[1],
                                  'gcpneg1_pval' = lcv$pval.fullycausal[2],
                                  'disease1_heritability' = lcv$h2.zscore[1],
                                  'disease2_heritability' = lcv$h2.zscore[2]))
}

# process
# filtering as done in original methods paper
lcv_results_processed <- lcv_results %>%
  filter(disease1_heritability > 7, disease2_heritability > 7, abs(corr_est) <= 0.9)
lcv_results_processed$gcp0_fdr <- p.adjust(lcv_results_processed$gcp0_pval, method = 'BH')
lcv_results_processed <- lcv_results_processed %>%
  mutate(assoc_label = case_when(gcp0_fdr < 0.01 ~ 'strictly_sig',  # FDR threshold used in methods paper
                                 gcp0_fdr >= 0.01 & gcp0_fdr < 0.05 ~ 'sig',
                                 .default = 'insig')) %>%
  mutate(causality_label = case_when(abs(gcp_est) > 0.6 ~ 'full_nearly_full_causality', .default = 'partial_causality'))

# export
write_csv(lcv_results, '~/data/internal/genomics/causality/lcv_results_all.csv')
write_csv(lcv_results_processed, '~/data/internal/genomics/causality/lcv_results_processed.csv')
