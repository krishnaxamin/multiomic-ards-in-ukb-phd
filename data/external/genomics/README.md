# Data present
- `chr_coordinates_hg19.csv`: chromosome lengths and coordinates on a continuous scale, on the hg19 assembly
- `berisa_pickrell_ld_blocks.csv`: co-ordinates (on hg19/GRch37) of pre-determined LD blocks, sourced from Berisa and Pickrell, 2016.

# Data absent
Participant-level data removed
- `array-genotyping_phenotype_file.csv`: phenotypes for the participants in the Genomics cohort constructed on the Cohort Browser.

Other data not present
- `ukb_snp_qc.txt`: from UK Biobank Resource 1955 (https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955).
- `ukb_geno_ref_alt_info.csv`: REF/ALT alleles of genotyped variants in the UKB, derived from `ukb_snp_qc.txt` using `~/ch2_genomics/get_genotyped_ref_alt.py`.
- `ukb_imputed_data/*_info_passed.mfi.txt`: information on imputed variants, covering (1) unique ID (2) snpID (e.g., rsID or Affymetrix ID) (3) position (4) REF/ALT (5) MAF (6) INFO score
- `UKBB.EUR.l2.ldscore`: LD scores for the Pan-UKB 'EUR' population, used in LCV modelling. Downloaded via https://pan-dev.ukbb.broadinstitute.org/docs/ld/index.html on 19 July, 2025. 