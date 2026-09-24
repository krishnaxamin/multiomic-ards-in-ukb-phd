"""
Use UK Biobank Resource 1955 (https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955) to get REF/ALT alleles for the
genotyped variants in the UK Biobank.
"""
from pandas import read_csv

info = read_csv('~/data/external/genomics/ukb_snp_qc.txt', delim_whitespace=True)

info_slimmed = info[['rs_id', 'chromosome', 'position', 'allele1_ref', 'allele2_alt']]

info_slimmed.columns = ['id', 'chr', 'pos', 'ref', 'alt']

info_slimmed.to_csv('~/data/external/genomics/ukb_geno_ref_alt_info.csv', index=False)
