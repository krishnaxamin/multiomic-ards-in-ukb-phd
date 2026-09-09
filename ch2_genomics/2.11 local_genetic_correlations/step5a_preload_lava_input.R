# STEP 5A: pre-load LAVA input object (avoids continual loading for each array task)

library(LAVA)

input <- process.input(input.info.file = '~/ch2_genomics/2.11 local_genetic_correlations/int_data/input.info.txt',           # input info file
                       sample.overlap.file = '~/ch2_genomics/2.11 local_genetic_correlations/int_data/sample_overlap.txt',   # sample overlap file (can be set to NULL if there is no overlap)
                       ref.prefix = '~/data/external/genomics/lava/ld_ref/lava-ukb-v1.1',                    # reference data prefix
                       phenos = NULL)  # NULL = all phenotypes to be processed

saveRDS(input, '~/ch2_genomics/2.11 local_genetic_correlations/int_data/lava_input.rds')
