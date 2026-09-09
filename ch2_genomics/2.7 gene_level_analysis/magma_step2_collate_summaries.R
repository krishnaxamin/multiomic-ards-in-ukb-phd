# Collate GWAS summary statistics for a disease and format for MAGMA.
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

disease_field <- commandArgs(trailingOnly = TRUE)[1]

adjustment <- commandArgs(trailingOnly = TRUE)[3]
if (adjustment == 'minimal'){
  svat_suffix = ''
  magma_suffix = 'minimal'
} else if (adjustment == 'lifestyles') {
  svat_suffix = '_lifestyles'
  magma_suffix = 'lifestyles'
}

print(paste0('Disease: ', disease_field, '. Adjustment: ', adjustment))
flush.console()
Sys.sleep(1)

# combine summary stats for all variants in disease
disease_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat', svat_suffix, '/', disease_field, '/', disease_field, '_pan-ukbb-eur.geno.svat'),
                         show_col_types = FALSE, num_threads = 1)
disease_imputed_files <- list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat', svat_suffix, '/', disease_field, '/'),
                                    pattern = 'imputed.svat', full.names = TRUE)
disease_imputed <- data.frame()
for (file in disease_imputed_files){
  if (!grepl('imputed.svat.index', file)){
    disease_imputed <- rbind(disease_imputed, read_tsv(file, show_col_types = FALSE, num_threads = 1) %>% 
                               filter(imputationInfo > 0.8) %>%
                               filter(abs(AF_Allele2) - 0.5 < 0.49))
  }
}

# for MAGMA, only need SNP ID and pval
disease_data <- rbind(disease_geno %>%
                        select(MarkerID, CHR, POS, Allele2, Allele1, p.value) %>%
                        rename(snpid = MarkerID, chr = CHR, bpos = POS, a2 = Allele1, a1 = Allele2, pval = p.value),  
                      disease_imputed %>%
                        select(MarkerID, CHR, POS, Allele2, Allele1, p.value) %>%
                        rename(snpid = MarkerID, chr = CHR, bpos = POS, a2 = Allele1, a1 = Allele2, pval = p.value)) %>%
  rowwise() %>%
  mutate(alleles_sorted = paste(sort(c(a1, a2)), collapse = '_')) %>% 
  ungroup() %>%
  mutate(uniqueid_sorted_alleles = paste(chr, bpos, alleles_sorted, sep = '_')) %>%  # generate uniqueID with sorted allele order
  distinct(uniqueid_sorted_alleles, .keep_all = TRUE) %>%  # sort on this to account for genotyped and imputed versions of the same variant having different allele orders
  # select(-alleles_sorted, -uniqueid_sorted_alleles) %>%
  select(snpid, pval)
  # distinct(snpid, a1, a2, .keep_all = TRUE)

write_tsv(disease_data, paste0('~/ch2_genomics/2.7 gene_level_analysis/gwas_summary_stats_for_magma/', disease_field, '_magma_', magma_suffix, '.tsv'))
