library(devtools)
# install_github("kpmainali/CooccurrenceAffinity")
library(CooccurrenceAffinity)
library(tidyverse)

# Calculate co-occurrence measure for disease pairs in a given cohort

# Genomics ------

# get number of cases per disease pair
pairwise_counts <- read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences.csv')

# get number of cases per disease
num_ppl_per_disease_df <- read_csv('~/data/internal/genomics/genomics_qc2_pan_ukbb_eur_common-unisex-disease_num-diagnosed_case-control.csv')

# two vectors, with each vector elementwise containing one half of a disease pair
disease1_names <- pairwise_counts$disease1
disease2_names <- pairwise_counts$disease2

num_indivs <- 423223 # N in CooccurrenceAffinity::ML.alpha()

alpha_hats <- rep(0, length(disease1_names))
alpha_hat_cis_lower <- rep(0, length(disease1_names))
alpha_hat_cis_upper <- rep(0, length(disease1_names))
alpha_hat_pvals <- rep(0, length(disease1_names))
for (i in 1:length(disease1_names)){
  pair1 <- disease1_names[i]
  pair2 <- disease2_names[i]
  pair1_pair2_ill <- as.numeric(pairwise_counts[i, 'cooccurrence_count']) # X in CooccurrenceAffinity::ML.alpha()
  pair1_ill_total <- num_ppl_per_disease_df[num_ppl_per_disease_df$disease_field == pair1, ]$num_cases # mA in CooccurrenceAffinity::ML.alpha()
  pair2_ill_total <- num_ppl_per_disease_df[num_ppl_per_disease_df$disease_field == pair2, ]$num_cases # mB in CooccurrenceAffinity::ML.alpha()
  cooccurrence_info <- ML.Alpha(pair1_pair2_ill, c(pair1_ill_total, pair2_ill_total, num_indivs), lev = 0.95)
  alpha_hats[i] <- cooccurrence_info$est # max. likelihood estimate of alpha
  alpha_hat_cis_lower[i] <- cooccurrence_info$CI.Blaker[1] # lower bound for 95% C.I. (Blaker)
  alpha_hat_cis_upper[i] <- cooccurrence_info$CI.Blaker[2] # upper bound for 95% C.I. (Blaker)
  alpha_hat_pvals[i] <- cooccurrence_info$pval # p-value for testing H0: alpha=0. By Blaker method (for compatibility with Blaker C.I.)
}

pairwise_counts$alpha_hat <- alpha_hats
pairwise_counts$alpha_hat_ci_lower <- alpha_hat_cis_lower
pairwise_counts$alpha_hat_ci_upper <- alpha_hat_cis_upper
pairwise_counts$alpha_hat_pval <- alpha_hat_pvals

pairwise_counts$alpha_hat_abs <- abs(pairwise_counts$alpha_hat)

write.csv(pairwise_counts, '~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv', row.names = FALSE)

# Omics Function ------
calculate_mle_omics <- function(pairwise_counts, phenotypes_df){
  
  # get number of cases per disease
  num_ppl_per_disease_df <- data.frame(disease_field = colnames(phenotypes_df %>%
                                                                  select(-eid)),
                                       num_cases = colSums(phenotypes_df %>%
                                                             select(-eid))) %>%
    as.data.frame(row.names = 1:nrow(.))
  
  # two vectors, with each vector elementwise containing one half of a disease pair
  disease1_names <- pairwise_counts$disease1
  disease2_names <- pairwise_counts$disease2
  
  num_indivs <- nrow(phenotypes_df) # N in CooccurrenceAffinity::ML.alpha()
  
  # CHECK FROM HERE ONWARDS
  alpha_hats <- rep(0, length(disease1_names))
  alpha_hat_cis_lower <- rep(0, length(disease1_names))
  alpha_hat_cis_upper <- rep(0, length(disease1_names))
  alpha_hat_pvals <- rep(0, length(disease1_names))
  for (i in 1:length(disease1_names)){
    pair1 <- disease1_names[i]
    pair2 <- disease2_names[i]
    pair1_pair2_ill <- as.numeric(pairwise_counts[i, 'cooccurrence_count']) # X in CooccurrenceAffinity::ML.alpha()
    pair1_ill_total <- num_ppl_per_disease_df[num_ppl_per_disease_df$disease_field == pair1, ]$num_cases # mA in CooccurrenceAffinity::ML.alpha()
    pair2_ill_total <- num_ppl_per_disease_df[num_ppl_per_disease_df$disease_field == pair2, ]$num_cases # mB in CooccurrenceAffinity::ML.alpha()
    if (pair1_pair2_ill * pair1_ill_total * pair2_ill_total == 0){
      # if any of the counts is 0, then calculating co-occurrence is irrelevant
      alpha_hats[i] <- -100 
      alpha_hat_cis_lower[i] <- -100 
      alpha_hat_cis_upper[i] <- -100
      alpha_hat_pvals[i] <- 1 
      next
    }
    cooccurrence_info <- ML.Alpha(pair1_pair2_ill, c(pair1_ill_total, pair2_ill_total, num_indivs), lev = 0.95)
    alpha_hats[i] <- cooccurrence_info$est # max. likelihood estimate of alpha
    alpha_hat_cis_lower[i] <- cooccurrence_info$CI.Blaker[1] # lower bound for 95% C.I. (Blaker)
    alpha_hat_cis_upper[i] <- cooccurrence_info$CI.Blaker[2] # upper bound for 95% C.I. (Blaker)
    alpha_hat_pvals[i] <- cooccurrence_info$pval # p-value for testing H0: alpha=0. By Blaker method (for compatibility with Blaker C.I.)
  }
  
  pairwise_counts$alpha_hat <- alpha_hats
  pairwise_counts$alpha_hat_ci_lower <- alpha_hat_cis_lower
  pairwise_counts$alpha_hat_ci_upper <- alpha_hat_cis_upper
  pairwise_counts$alpha_hat_pval <- alpha_hat_pvals
  
  pairwise_counts$alpha_hat_abs <- abs(pairwise_counts$alpha_hat)
  
  return(pairwise_counts)
  
}

# Proteomics, phenotypes prior to blood taking ------

# get number of cases per disease pair
pairwise_counts <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences.csv', show_col_types = FALSE)

# get phenotypes
prior_phenotypes_df <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv', show_col_types = FALSE)

# run function
result <- calculate_mle_omics(pairwise_counts, phenotypes_df)
write_csv(result, '~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences_alphamle.csv')

# Proteomics, full cohort ------

# get number of cases per disease pair
pairwise_counts <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences.csv', show_col_types = FALSE)

# get phenotypes
phenotypes_df <- read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_phenotypes.csv', show_col_types = FALSE)

# run function
result <- calculate_mle_omics(pairwise_counts, phenotypes_df)
write_csv(result, '~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
