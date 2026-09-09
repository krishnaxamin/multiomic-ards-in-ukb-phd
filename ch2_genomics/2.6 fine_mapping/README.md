#  Code
- `2.6.1 step1_identify_blocks.py`: use results from minimally adjusted and lifestyle-adjusted GWAS to construct LD-independent blocks for (1) LD statistic calculations (2) fine mapping (3) colocalisation. Corresponds to **Methods M2.6.1**.
- `2.6.2 step2_ld_calcs.sh`: Calculate LD statistics (genotypic correlation, _r_) for the LD-independent blocks in `./blocks_for_ld_calcs.csv`. Corresponds to **Methods M2.6.2**.
- `2.6.3 fine_mapping/`: Directory containing scripts to perform fine mapping and analyse results
- `plot_fine_mapping_dists_counts.py`: Generate plots and datasets (used to make tables) for fine mapping results, for use in the thesis.

# Data present
- `blocks_for_ld_calcs.csv`: LD-independent blocks, where each block contains variants associated with ARDs through both minimally adjusted and lifestyle-adjusted GWAS, and neighbouring variants within LD (taken as within 1Mb). LD statistics (signed genotypic correlation, _r_) calculated within these blocks.
- `disease_ld_indep_blocks_minimal.csv` & `disease_ld_indep_blocks_lifestyles.csv`: LD-independent blocks, where each block contains variants associated with a single ARD through one of minimally adjusted or lifestyle-adjusted GWAS. Used in fine mapping.