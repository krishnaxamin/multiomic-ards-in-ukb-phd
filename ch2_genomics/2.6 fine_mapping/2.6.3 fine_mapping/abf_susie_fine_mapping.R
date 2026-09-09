# fine mapping using `coloc` general script
# run ABF fine mapping for all diseases for which fine mapping not yet complete in the coloc pipeline
# export results
# similarly run SuSiE fine mapping if possible, export if so

if (!require("coloc")) install.packages("coloc")
library(coloc)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


region_index <- commandArgs(trailingOnly = TRUE)[1]
adjustment <- commandArgs(trailingOnly = TRUE)[2]  # minimal OR lifestyles

svat_suffix = ''
if (adjustment == 'lifestyles') { svat_suffix = '_lifestyles' }

disease_regions <- read.csv(paste0('~/ch2_genomics/2.6 fine_mapping/2.6.3 fine_mapping/disease_ld_indep_blocks_', adjustment, '.csv'))
regions_lded <- read.csv('~/ch2_genomics/2.6 fine_mapping/2.6.3 fine_mapping/blocks_for_ld_calcs.csv')

region <- disease_regions[region_index, ]

chr <- region$chr
start <- region$start
end <- region$end
disease <- region$disease

# check whether fine mapping done already in coloc pipeline
if (file.exists(paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/abf/robjs/', disease, '_', chr, '_', start, '_', end, '_abf.rds'))) {
  stop(paste0(disease, ' for region ', chr, ':', start, '-', end, ' has already been fine mapped.'))
}

disease_info <- read.csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_field <- (disease_info %>%
                    filter(icd10_three_letter == disease))$disease_field

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

# read in disease results data
disease_geno <- read_tsv(paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat', svat_suffix, '/', disease_field, '/', disease_field, '_pan-ukbb-eur.geno.svat'),
                         show_col_types = FALSE, num_threads = 1)
disease_imputed_files <- as.list(list.files(path = paste0('~/ch2_genomics/2.5 full_gwas_runs/step2_svat', svat_suffix, '/', disease_field, '/'),
                                            pattern = 'imputed.svat', full.names = TRUE))
disease_imputed <- data.frame()
for(i in seq_along(l=disease_imputed_files)){
  if (!(grepl('index', disease_imputed_files[[i]]))){
    disease_imputed <- rbind(disease_imputed, 
                             read_tsv(disease_imputed_files[[i]], 
                                      show_col_types = FALSE, num_threads = 1) %>% filter(imputationInfo > 0.8))
  }
}
disease_data <- rbind(disease_geno %>%
                        select(CHR, POS, BETA, SE, Allele1, Allele2),
                      disease_imputed %>%
                        select(CHR, POS, BETA, SE, Allele1, Allele2)) %>%
  mutate(uniqueID = paste0(CHR, ':', POS, ':', Allele1, ':', Allele2)) %>%
  mutate(var = SE ^ 2) %>%
  distinct(uniqueID, .keep_all = TRUE) %>%
  filter(CHR == chr, POS >= start, POS <= end) %>%
  filter(uniqueID %in% variants_to_include) %>%
  mutate(order_col = match(uniqueID, variants_to_include)) %>%
  arrange(order_col) %>%
  select(-order_col, -SE)

# set up coloc objects

# set up LD matrices. These can still have NaNs in them, so iteratively remove rows and columns with the most NaNs until no NaNs remain 
disease_ld <- ld_data_cleaned[disease_data$uniqueID, disease_data$uniqueID]
num_nans <- unique(rowSums(is.na(disease_ld)))
while (length(num_nans) > 1){
  disease_ld <- disease_ld[rowSums(is.na(disease_ld)) < max(num_nans), colSums(is.na(disease_ld)) < max(num_nans)]
  num_nans <- unique(rowSums(is.na(disease_ld)))
}

vars_non_nan_ld <- colnames(disease_ld)

# update disease_data to reflect new LD matrices
disease_data_updated <- disease_data %>% filter(uniqueID %in% vars_non_nan_ld)

disease_coloc <- list(beta = disease_data_updated$BETA,
                      varbeta = disease_data_updated$var,
                      N = 423223,
                      type = 'cc',
                      s = 1/((disease_info %>%
                                filter(icd10_three_letter == {{disease}}))$case_control_ratio_genomics + 1),
                      LD = as.matrix(disease_ld),
                      snp = disease_data_updated$uniqueID,
                      position = disease_data_updated$POS)

# ABF on disease
disease_abf <- finemap.abf(disease_coloc)
saveRDS(disease_abf, paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/abf/robjs/', disease, '_', chr, '_', start, '_', end, '_abf.rds'))

# try SuSiE
tryCatch(
  expr = {
    disease_susie <- runsusie(disease_coloc)
    saveRDS(disease_susie, paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/susie/robjs/', disease, '_', chr, '_', start, '_', end, '_susie.rds'))
  },
  error = function(e) {
    print(paste0('For ', disease, ':'))
    print(e)
  }
)
