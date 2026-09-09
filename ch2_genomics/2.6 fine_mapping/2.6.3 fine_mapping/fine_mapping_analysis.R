# Collate results from fine mapping using SuSiE or "normal" fine mapping and identify causal variants for each disease of varying degrees of confidence
# Confidence is determined by the p(causal) value (PP or pip)

library(tidyverse)

adjustment <- 'lifestyles'

# get paths to SuSiE and ABF objects
susie_objects <- list.files(path = paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/susie/robjs/'),
                            pattern = 'susie.rds', full.names = TRUE)
abf_objects <- list.files(path = paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/abf/robjs/'),
                          pattern = 'abf.rds', full.names = TRUE)

# set up stores for variants
# peak = most credible variant in each credible set from SuSiE fine mappings & top variant from ABF fine mappings
# all = all credible/fine mapped variants
peak_fine_mapped_variants <- data.frame()
all_fine_mapped_variants <- data.frame()

# pull from SuSiE fine mappings
for (susie_object in susie_objects){
  susie <- readRDS(susie_object)
  
  # get disease and other info
  info_from_file_name <- str_split_1(str_split_1(susie_object, '/')[8], '_')
  disease <- info_from_file_name[1]
  chr <- info_from_file_name[2]
  region_start <- info_from_file_name[3]
  region_end <- info_from_file_name[4]
  
  
  # access credible sets
  susie_credible_sets <- susie$sets$cs
  
  # want to pull out the most credible variant from each credible set
  credible_set_counter <- 0
  for (cs in susie_credible_sets){
    credible_set_counter <- credible_set_counter + 1
    cs_pips <- susie$pip[cs]
    
    most_credible_variant <- data.frame(disease = disease, 
                                        variant = names(sort(cs_pips, decreasing = TRUE))[1],
                                        pip = unname(sort(cs_pips, decreasing = TRUE))[1],
                                        algo = 'SuSiE',
                                        chr = chr,
                                        ld_region_start = region_start,
                                        ld_region_end = region_end,
                                        disease_credible_set = credible_set_counter)
    peak_fine_mapped_variants <- rbind(peak_fine_mapped_variants, most_credible_variant)
    
    credible_variants <- data.frame(disease = disease,
                                    variant = names(sort(cs_pips, decreasing = TRUE)),
                                    pip = unname(sort(cs_pips, decreasing = TRUE)),
                                    algo = 'SuSiE',
                                    chr = chr,
                                    ld_region_start = region_start,
                                    ld_region_end = region_end,
                                    disease_credible_set = credible_set_counter)
    all_fine_mapped_variants <- rbind(all_fine_mapped_variants, credible_variants)
  }
}

# pull from ABF fine mappings
for (abf_object in abf_objects){
  abf <- readRDS(abf_object)
  
  # get disease and other info
  info_from_file_name <- str_split_1(str_split_1(abf_object, '/')[8], '_')
  disease <- info_from_file_name[1]
  chr <- info_from_file_name[2]
  region_start <- info_from_file_name[3]
  region_end <- info_from_file_name[4]
  
  # only incorporate these fine mapping results if the disease was not 
  # if ((disease %in% non_susied_diseases) & (region_start %in% non_susied_regions$start) & (region_end %in% non_susied_regions$end)){
  abf_variants <- abf %>%
    filter(!is.na(position)) %>%
    select(snp, SNP.PP) %>%
    rename(variant = snp, pip = SNP.PP) %>%
    mutate(disease = {{disease}},
           algo = 'ABF',
           chr = {{chr}},
           ld_region_start = {{region_start}},
           ld_region_end = {{region_end}},
           disease_credible_set = 0)
  all_fine_mapped_variants <- rbind(all_fine_mapped_variants, abf_variants)
  
  abf_top_variant <- (abf_variants %>%
                        arrange(desc(pip)))[1, ]
  peak_fine_mapped_variants <- rbind(peak_fine_mapped_variants, abf_top_variant)
  # } 
}

# export fine mapped variants
write_csv(all_fine_mapped_variants, paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/all_fine_mapped_variants.csv'))
write_csv(peak_fine_mapped_variants, paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/peak_fine_mapped_variants.csv'))

# identify high/medium-confidence top causal variants
high_conf_peak_causal <- peak_fine_mapped_variants %>% filter(pip > 0.9)
med_conf_peak_causal <- peak_fine_mapped_variants %>% filter(pip > 0.5)

# identify high/medium-confidence causal variants from whole set
high_conf_causal <- all_fine_mapped_variants %>% filter(pip > 0.9)
med_conf_causal <- all_fine_mapped_variants %>% filter(pip > 0.5)

# export peak causal variants with pip > 0.5: lenient threshold
write_csv(med_conf_peak_causal, paste0('~/data/internal/genomics/fine_mapping/', adjustment, '/peak_pip0.5_fine_mapped_variants.csv'))
