# STEP 5B: run LAVA for genes that are in LD blocks for which there is FDR permissive < 0.05 correlation

# if (!require("remotes", quietly=T)) install.packages("remotes")
# remotes::install_github("josefin-werme/LAVA")
library(LAVA)
library(tidyverse)

# Read in data ------
locus_idx <- commandArgs(trailingOnly = TRUE)[1]
## Read in summary statistics and related info ------

input <- readRDS('~/ch2_genomics/2.11 local_genetic_correlations/int_data/lava_input.rds')

## Read in locus info file ------
loci <- read.loci('~/ch2_genomics/2.11 local_genetic_correlations/int_data/loci_follow_up_genes.txt')
locus_name <- loci[locus_idx, 'LOC']
n_loci <- nrow(loci)
# head(loci)  # inspect the locus file

print(paste0(Sys.time(), ': Analysing ', locus_name, ' (', locus_idx, '/', n_loci, ')'))
flush.console()
Sys.sleep(1)

# Create a locus object for the ith locus to prepare it for analysis
locus = process.locus(loci[locus_idx, ], input)

# ls(locus)                               # inspect locus
# c(locus$chr, locus$start, locus$stop)   # locus coordinates
# str(locus$snps)                         # locus snps
# locus$n.snps                            # N snps
# locus$omega                             # genetic covariance matrix
# locus$omega.cor                         # standardised genetic covariance matrix
# locus$phenos                            # locus phenotypes

# Testing ------
# p < 0.05 may allow false positives, but these will likely lead to null bivariate tests, 
# so better to be generous and power-aware, and multiple testing-correct at the bivariate test
locus_result <- run.univ.bivar(locus, univ.thresh = 0.05)  

# Export ------
write_csv(locus_result$univ %>%
            select(phen, h2.obs, p) %>%
            mutate(locus = locus_name) %>%
            rename(disease = phen, pval = p),
          paste0('~/ch2_genomics/2.11 local_genetic_correlations/results/heritability/', locus_name, '_lava_heritability.csv'))

if (is.null(locus_result$bivar)) { 
  print(paste0(Sys.time(), ': no bivariate result.'))
  flush.console()
  Sys.sleep(1)
} else {
  write_csv(locus_result$bivar %>%
              select(!contains('r2')) %>%
              mutate(locus = locus_name) %>%
              rename(disease1 = phen1, disease2 = phen2, pval = p),
            paste0('~/ch2_genomics/2.11 local_genetic_correlations/results/correlation/', locus_name, '_lava_correlations.csv'))
}

