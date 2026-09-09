# coloc general script
# run ABF fine mapping and coloc for all disease pairs
# export results
# run SuSiE fine mapping and coloc to see if possible, export if so

# from author of coloc, who falls back on ABF if SuSiE doesn't run
#  https://github.com/chr1swallace/coloc/issues/90

# Performs SuSiE fine mapping and colocalisation for disease pairs that share associated variants. 
# For each disease pair, performs fine mapping and colocalisation for each region (defined in steps 1 and 2).
# Each region contains within them variants associated with either disease or both diseases that are <1Mb from the nearest relevant associated variant.

if (!require("coloc")) install.packages("coloc")
library(coloc)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

# read in data ------

region_index <- commandArgs(trailingOnly = TRUE)[1]

combo_regions <- read.csv(paste0('~/ch2_genomics/2.13 colocalisations/combo_regions_lifestyles.csv'))
regions_lded <- read.csv('~/ch2_genomics/2.6 fine_mapping/blocks_for_ld_calcs.csv')

region <- combo_regions[region_index, ]

chr <- region$chr
start <- region$start
end <- region$end
disease1 <- str_split(region$combo, '-')[[1]][1]
disease2 <- str_split(region$combo, '-')[[1]][2]

disease_info <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease1_field <- (disease_info %>%
                     filter(icd10_three_letter == disease1))$disease_field
disease2_field <- (disease_info %>%
                     filter(icd10_three_letter == disease2))$disease_field

# determine which LD region contains the desired region
ld_region <- regions_lded %>%
  filter(chr == {{chr}}, start <= {{start}}, end >= {{end}})

# read in LD data 
ld_variant_ids <- read_tsv(paste0('~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_', ld_region$chr, '_', ld_region$start, '_', ld_region$end, '.snplist'),
                           show_col_types = FALSE, col_names = FALSE, num_threads = 1)$X1
ld_data <- read_tsv(paste0('~/data/internal/genomics/ld_calcs/ld_mats/ld_mat_', ld_region$chr, '_', ld_region$start, '_', ld_region$end, '.ld'),
                    show_col_types = FALSE, col_names = ld_variant_ids, num_threads = 1)
ld_data <- data.frame(ld_data)
rownames(ld_data) <- ld_variant_ids
colnames(ld_data) <- ld_variant_ids

# remove rows and columns with all NANs
ld_data_cleaned <- ld_data[rowSums(is.na(ld_data)) < ncol(ld_data), colSums(is.na(ld_data)) < nrow(ld_data)]

# get variants which have LD data
variants_to_include <- rownames(ld_data_cleaned)

# read in disease1 results data
disease1_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease1_field, '/', disease1_field, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1)
disease1_imputed_files <- as.list(list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease1_field, '/'),
                                             pattern = 'imputed.svat', full.names = TRUE))
disease1_imputed <- data.frame()
for(i in seq_along(l=disease1_imputed_files)){
  if (!(grepl('index', disease1_imputed_files[[i]]))){
    disease1_imputed <- rbind(disease1_imputed, 
                              read_tsv(disease1_imputed_files[[i]], 
                                       show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease1_data <- rbind(disease1_geno %>%
                         select(CHR, POS, BETA, SE, Allele1, Allele2),
                       disease1_imputed %>%
                         select(CHR, POS, BETA, SE, Allele1, Allele2)) %>%
  mutate(uniqueID = paste0(CHR, ':', POS, ':', Allele1, ':', Allele2)) %>%
  mutate(var = SE ^ 2) %>%
  distinct(uniqueID, .keep_all = TRUE) %>%
  filter(CHR == chr, POS >= start, POS <= end) %>%
  filter(uniqueID %in% variants_to_include) %>%
  mutate(order_col = match(uniqueID, variants_to_include)) %>%  # match(x, y) finds position of x in y
  arrange(order_col) %>%
  select(-order_col, -SE)

# read in disease2 results data
disease2_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/', disease2_field, '_pan-ukbb-eur.geno.svat'),
                          show_col_types = FALSE, num_threads = 1)
disease2_imputed_files <- as.list(list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/', disease2_field, '/'),
                                             pattern = 'imputed.svat', full.names = TRUE))
disease2_imputed <- data.frame()
for(i in seq_along(l=disease2_imputed_files)){
  if (!(grepl('index', disease2_imputed_files[[i]]))){
    disease2_imputed <- rbind(disease2_imputed, 
                              read_tsv(disease2_imputed_files[[i]], 
                                       show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease2_data <- rbind(disease2_geno %>%
                         select(CHR, POS, BETA, SE, Allele1, Allele2),
                       disease2_imputed %>%
                         select(CHR, POS, BETA, SE, Allele1, Allele2)) %>%
  mutate(uniqueID = paste0(CHR, ':', POS, ':', Allele1, ':', Allele2)) %>%
  mutate(var = SE ^ 2) %>%
  distinct(uniqueID, .keep_all = TRUE) %>%
  filter(CHR == chr, POS >= start, POS <= end) %>%
  filter(uniqueID %in% variants_to_include) %>%
  mutate(order_col = match(uniqueID, variants_to_include)) %>%
  arrange(order_col) %>%
  select(-order_col, -SE)

if (nrow(disease1_data) != nrow(disease2_data)){
  print('disease1_data and disease2_data do not have the same number of variants')
}  # this should be TRUE

# set up coloc objects ------

# set up LD matrices. These can still have NaNs in them, so iteratively remove rows and columns with the most NaNs until no NaNs remain 
disease1_ld <- ld_data_cleaned[disease1_data$uniqueID, disease1_data$uniqueID]
num_nans <- unique(rowSums(is.na(disease1_ld)))
while (length(num_nans) > 1){
  disease1_ld <- disease1_ld[rowSums(is.na(disease1_ld)) < max(num_nans), colSums(is.na(disease1_ld)) < max(num_nans)]
  num_nans <- unique(rowSums(is.na(disease1_ld)))
}

vars_non_nan_ld <- colnames(disease1_ld)

disease2_ld <- ld_data_cleaned[vars_non_nan_ld, vars_non_nan_ld]

# update disease1/2_data to reflect new LD matrices
disease1_data_updated <- disease1_data %>% filter(uniqueID %in% vars_non_nan_ld)
disease2_data_updated <- disease2_data %>% filter(uniqueID %in% vars_non_nan_ld)

disease1_coloc <- list(beta = disease1_data_updated$BETA,
                       varbeta = disease1_data_updated$var,
                       N = 423223,
                       type = 'cc',
                       s = 1/((disease_info %>%
                                 filter(icd10_three_letter == {{disease1}}))$case_control_ratio_genomics + 1),
                       LD = as.matrix(disease1_ld),
                       snp = disease1_data_updated$uniqueID,
                       position = disease1_data_updated$POS)
disease2_coloc <- list(beta = disease2_data_updated$BETA,
                       varbeta = disease2_data_updated$var,
                       N = 423223,
                       type = 'cc',
                       s = 1/((disease_info %>%
                                 filter(icd10_three_letter == {{disease2}}))$case_control_ratio_genomics + 1),
                       LD = as.matrix(disease2_ld),
                       snp = disease2_data_updated$uniqueID,
                       position = disease2_data_updated$POS)

# abf ------

coloc <- coloc.abf(disease1_coloc, disease2_coloc)
saveRDS(coloc, paste0('~/ch2_genomics/2.13 colocalisations/results/abf/robjs/',  disease1, '_', disease2, '_', chr, '_', start, '_', end, '_abf_coloc.rds'))

print('ABF fine mapping and coloc complete.')

## export abf ------

# export all coloc results, irrespective of PP.H4
coloc_summary <- data.frame(as.list(coloc$summary))
if (nrow(coloc_summary) > 0){
  coloc_summary$combo <- region$combo
  coloc_summary$region <- paste(region$chr, region$start, region$end, sep = '-')
  write_csv(coloc_summary, paste0('~/ch2_genomics/2.13 colocalisations/results/abf/results/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_abf_coloc_summary.csv'))
}

coloc_results <- coloc$results
if (nrow(coloc_results) > 0){
  coloc_results$combo <- region$combo
  coloc_results$region <- paste(region$chr, region$start, region$end, sep = '-')
  write_csv(coloc_results, paste0('~/ch2_genomics/2.13 colocalisations/results/abf/results/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_abf_coloc_results.csv'))
}

# plot sensitivity for row which passes the PP.H4 filter
coloc_hits <- coloc_results %>% filter(SNP.PP.H4 > 0.9)
print(paste0(nrow(coloc_hits), ' variants pass SNP.PP.H4 > 0.9'))
if (nrow(coloc_hits) > 0){
  png(paste0('~/ch2_genomics/2.13 colocalisations/results/abf/sensitivity_plots/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_abf_coloc_sensitivity.png'),
      width = 1123, height = 794, res = 96)
  sensitivity(coloc,"H4 > 0.9")
  dev.off()
  pdf(paste0('~/ch2_genomics/2.13 colocalisations/results/abf/sensitivity_plots/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_abf_coloc_sensitivity.pdf'),
      width = 11.693, height = 8.268)
  sensitivity(coloc,"H4 > 0.9")
  dev.off()
}

# susie ------

susie_success_counter = 0
# try SuSiE on disease1
tryCatch(
  expr = {
    disease1_susie <- runsusie(disease1_coloc)
    susie_success_counter = susie_success_counter + 1
    print(paste0('SuSiE on ', disease1, ' successful.'))
    },
  error = function(e) {
    print(paste0('For ', disease1, ':'))
    print(e)
  }
)

# try SuSiE on disease2
tryCatch(
  expr = {
    disease2_susie <- runsusie(disease2_coloc)
    susie_success_counter = susie_success_counter + 1
    print(paste0('SuSiE on ', disease2, ' successful.'))
  },
  error = function(e) {
    print(paste0('For ', disease2, ':'))
    print(e)
  }
)

if (susie_success_counter < 2) {
  stop('SuSiE fine mapping was not completed successfully for both diseases, so SuSiE colocalisation cannot proceed.')
}

if (is.null(disease1_susie$sets$cs) | is.null(disease2_susie$sets$cs)){
  stop('SuSiE credible sets were not present for both diseases, so SuSiE coloc could not be run.')
}

# try SuSiE coloc 
coloc <- coloc.susie(disease1_susie, disease2_susie)
saveRDS(coloc, paste0('~/ch2_genomics/2.13 colocalisations/results/susie/robjs/',  disease1, '_', disease2, '_', chr, '_', start, '_', end, '_susie_coloc.rds'))

## export susie ------

# export all coloc results, irrespective of PP.H4

coloc_summary <- coloc$summary
if (nrow(coloc_summary) > 0){
  coloc_summary$combo <- region$combo
  coloc_summary$region <- paste(region$chr, region$start, region$end, sep = '-')
  write_csv(coloc_summary, paste0('~/ch2_genomics/2.13 colocalisations/results/susie/results/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_susie_coloc_summary.csv'))
}

coloc_results <- coloc$results
if (nrow(coloc_results) > 0){
  coloc_results$combo <- region$combo
  coloc_results$region <- paste(region$chr, region$start, region$end, sep = '-')
  write_csv(coloc_results, paste0('~/ch2_genomics/2.13 colocalisations/results/susie/results/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_susie_coloc_results.csv'))
}

print('SuSiE coloc complete and results exported.')

# plot sensitivity for rows which pass the PP.H4 filter
coloc_hits <- coloc_summary %>% filter(PP.H4.abf > 0.9)
print(paste0(nrow(coloc_hits), ' variants pass SNP.PP.H4 > 0.9'))
if (nrow(coloc_hits) > 0){
  
  # sensitivity plotting
  rows_to_plot_sensitivity <- (coloc_summary %>% mutate(row_number = rownames(coloc_summary)) %>% filter(PP.H4.abf > 0.9))$row_number
  for (row in rows_to_plot_sensitivity){
    png(paste0('~/ch2_genomics/2.13 colocalisations/results/susie/sensitivity_plots/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_susie_coloc', row, '_sensitivity.png'),
        width = 1123, height = 794, res = 96)
    sensitivity(coloc,"H4 > 0.9", row = as.numeric(row), dataset1 = disease1_coloc, dataset2 = disease2_coloc)
    dev.off()
    pdf(paste0('~/ch2_genomics/2.13 colocalisations/results/susie/sensitivity_plots/', disease1, '_', disease2, '_', chr, '_', start, '_', end, '_susie_coloc', row, '_sensitivity.pdf'),
        width = 11.693, height = 8.268)
    sensitivity(coloc,"H4 > 0.9", row = as.numeric(row), dataset1 = disease1_coloc, dataset2 = disease2_coloc)
    dev.off()
  }
  
}
