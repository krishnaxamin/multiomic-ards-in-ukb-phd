# STEP 3: generate sample overlap file from LDSC results
library(tidyverse)

# collate diseaseX-diseaseY results ------
ldsc_results <- list.files(path = '~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/results/',
                           pattern = '.log', full.names = TRUE)

ldsc_covs <- data.frame()
for (ldsc_result in ldsc_results) {
  ldsc <- readLines(ldsc_result)
  genetic_covariance_intercept_line <- ldsc[which(ldsc == 'Genetic Covariance') + 4]
  if (grepl('Intercept: ', genetic_covariance_intercept_line)) {
    combo_cov <- as.numeric(str_split_1(genetic_covariance_intercept_line, ' ')[2])
  } else {
    stop(paste0(combo, ' was not read successfully.'))
  }
  disease_pair_split <- str_split_1(str_replace(basename(ldsc_result), '.log', ''), '-')
  ldsc_covs <- rbind(ldsc_covs, data.frame('disease1' = disease_pair_split[1],
                                           'disease2' = disease_pair_split[2],
                                           'genetic_covariance' = combo_cov))
}

# collate diseaseX-diseaseX results ------
ldsc_same_results <- list.files(path = '~/ch2_genomics/2.11 local_genetic_correlations/int_data/diseaseX_diseaseX_ldsc_results/',
                                pattern = '.log', full.names = TRUE)
ldsc_same_covs <- data.frame()
for (ldsc_same_result in ldsc_same_results) {
  ldsc <- readLines(ldsc_same_result)
  genetic_covariance_intercept_line <- ldsc[which(ldsc == 'Genetic Covariance') + 4]
  if (grepl('Intercept: ', genetic_covariance_intercept_line)) {
    combo_cov <- as.numeric(str_split_1(genetic_covariance_intercept_line, ' ')[2])
  } else {
    stop(paste0(combo, ' was not read successfully.'))
  }
  disease_pair_split <- str_split_1(str_replace(basename(ldsc_same_result), '.log', ''), '-')
  ldsc_same_covs <- rbind(ldsc_same_covs, data.frame('disease1' = disease_pair_split[1],
                                                     'disease2' = disease_pair_split[2],
                                                     'genetic_covariance' = combo_cov))
}

# organise into LAVA matrix ------
phenotypes <- sort(unique(c(ldsc_covs$disease1, ldsc_covs$disease2)))

cov_mat <- rbind(ldsc_covs,
                 ldsc_covs %>%
                   rename(disease1 = disease2, disease2 = disease1),
                 ldsc_same_covs) %>%
  pivot_wider(names_from  = disease2, values_from = genetic_covariance) %>%
  column_to_rownames("disease1") %>%
  as.matrix()

if (!identical(sort(rownames(cov_mat)), sort(colnames(cov_mat)))) {
  stop('Row and column names are not the same in the covariance matrix.')
}
if (!identical(sort(rownames(cov_mat)), phenotypes)) {
  stop('Row and column names are not the same as the phenotype set.')
}

cov_mat <- cov_mat[phenotypes, phenotypes]

mat = round(cov2cor(cov_mat), 5) # standardise
write.table(mat, '~/ch2_genomics/2.11 local_genetic_correlations/int_data/sample_overlap.txt', quote=F) # save
