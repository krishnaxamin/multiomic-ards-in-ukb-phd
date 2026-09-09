# STEP 0 (local): generate loci file
# LOC, CHR, START, STOP

library(tidyverse)

# LD blocks from Berisa and Pickrell (same as used in Donertas and LD block sharing)
ld_blocks <- read.table('~/data/external/genomics/berisa_pickrell_ld_blocks.bed',
                        header = TRUE) %>%
  rename(CHR = chr, START = start, STOP = stop) %>%
  mutate(LOC = paste0('block', row_number() - 1), CHR = str_remove(CHR, 'chr'))

# export
write_tsv(ld_blocks, '~/ch2_genomics/2.11 local_genetic_correlations/int_data/loci_ld_blocks_only.txt')
