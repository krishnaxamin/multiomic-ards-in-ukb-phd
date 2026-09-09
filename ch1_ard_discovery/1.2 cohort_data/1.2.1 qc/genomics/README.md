# Genomics QC pipeline
- Step 1: variant missingness
	- Step 2:  TSV of high-missingness variants
- Step 3: sample missingness
	- Step 4: TSV of high-missingness samples (not present as contains UKB participant-level information)
- Step 5: variant minor allele frequency (MAF)
	- Step 6: TSV of low-MAF variants
- Step 7: Hardy-Weinberg equilibrium (HWE)
	- Step 8: CSV of number of variants excluded by HWE filtering for each ARD

Odd-numbered steps are performed using bash. Even-numbered steps are the output from processing the results of the bash runs in Python.
	
Variants that are high-missingness, low-MAF, or both are listed in `~/data/internal/genomics/high_missingness_low_maf_variants_panukbb_eur.tsv`.