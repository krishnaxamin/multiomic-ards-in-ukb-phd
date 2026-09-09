# Plot absolute z-scores for regions on c19 surrounding variants associated with G30 and E11,
#  to see if extreme z-scores exist in G30 but not in E11. 
#  Extreme z-scores whose neighbours don't show commensurate signal and that don't work with the LD matrix 
#  (unsure why the two don't work together) are a possible reason why SuSiE fails
#  to run with G30's c19 region, citing 'The estimated prior variance is unreasonably large',
#  but why E11's c19 region, which overlaps, does work
# Each region contains within them variants associated disease that are <1Mb from the nearest relevant associated variant.

# Plots show a large incoherent peak (|z| > 40) in G30, but for E11, more coherent peaks with |z| < 8

if (!require("coloc")) install.packages("coloc")
library(coloc)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)


region_index <- commandArgs(trailingOnly = TRUE)[1] # 75 for G30, 28 for E11
adjustment <- commandArgs(trailingOnly = TRUE)[2]  # lifestyles

disease_regions <- read.csv(paste0('/data/home/bty207/gwas/fine_mapping/saige/disease_ld_indep_blocks_', adjustment, '.csv'))
regions_lded <- read.csv('/data/home/bty207/gwas/blocks_for_ld_calcs.csv')

region <- disease_regions[region_index, ]

chr <- region$chr
start <- region$start
end <- region$end
disease <- region$disease

disease_info <- read.csv('/data/home/bty207/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_field <- (disease_info %>%
                    filter(icd10_three_letter == disease))$disease_field

# determine which LD region contains the desired region
ld_region <- regions_lded %>%
  filter(chr == {{chr}}, start <= {{start}}, end >= {{end}})

# read in LD data 
ld_variant_ids <- read_tsv(paste0('/data/WHRI-Phenogenomics/krishna/ld_calcs/ld_mats/ld_mat_', ld_region$chr, '_', ld_region$start, '_', ld_region$end, '.snplist'),
                           show_col_types = FALSE, col_names = FALSE, num_threads = 1)$X1
ld_data <- read_tsv(paste0('/data/WHRI-Phenogenomics/krishna/ld_calcs/ld_mats/ld_mat_', ld_region$chr, '_', ld_region$start, '_', ld_region$end, '.ld'),
                    show_col_types = FALSE, col_names = ld_variant_ids, num_threads = 1)
ld_data <- data.frame(ld_data)
rownames(ld_data) <- ld_variant_ids
colnames(ld_data) <- ld_variant_ids

# remove rows and columns with all NANs
ld_data_cleaned <- ld_data[rowSums(is.na(ld_data)) < ncol(ld_data), colSums(is.na(ld_data)) < nrow(ld_data)]

# get variants which have LD data
variants_to_include <- rownames(ld_data_cleaned)

# read in disease results data
disease_geno <- read_tsv(paste0('/data/WHRI-Phenogenomics/krishna/ukb_saige_sparse_svat_', adjustment, '/', disease_field, '/', disease_field, '_pan-ukbb-eur.geno.svat'),
                         show_col_types = FALSE, num_threads = 1)
disease_imputed_files <- as.list(list.files(path = paste0('/data/WHRI-Phenogenomics/krishna/ukb_saige_sparse_svat_', adjustment, '/', disease_field, '/'), 
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
  mutate(order_col = match(uniqueID, variants_to_include)) %>%  # match(x, y) finds position of x in y
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

# add z-scores and plot for all variants for this region and disease
disease_data_updated <- disease_data_updated %>% mutate(z = abs(BETA / sqrt(var)))
ggplot(disease_data_updated, aes(x = POS, y = z)) + 
  geom_point() + 
  labs(title = paste0(
    'Absolute z-score of variants from SAIGE GWAS of ', disease, ' in ', ld_region$chr, ':', 
    ld_region$start, '-', ld_region$end), y = '|z-score|')