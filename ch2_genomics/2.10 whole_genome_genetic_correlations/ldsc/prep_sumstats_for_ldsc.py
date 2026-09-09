""" Script to collate and process SAIGE results for LDSC munging """

import argparse
import os
import pandas as pd

""" Digest command line arguments """
parser = argparse.ArgumentParser(description='Reading in disease whose SAIGE results to prep.')
parser.add_argument('--disease', type=str, required=True, help='disease field')
parser.add_argument('--adjustment', type=str, required=True, help='lifestyles or minimal')
args = parser.parse_args()  # parse arguments

disease = args.disease

""" Read in genotype data """
geno = pd.read_csv(
    f"~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/{disease}/{disease}_pan-ukbb-eur.geno.svat",
    sep='\s+', usecols=['MarkerID', 'Allele2', 'Allele1', 'BETA', 'N_case', 'N_ctrl', 'AF_Allele2', 'p.value']
).assign(imputationInfo=1.0)

""" Read in imputed data """
imputed_files = [x.name for x in os.scandir(
    f"~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/{disease}"
) if x.is_file() and 'imputed.svat' in x.name and 'svat.index' not in x.name]

imputed = pd.DataFrame()
for file in imputed_files:
    imputed = pd.concat([imputed,
                         pd.read_csv(f"~/ch2_genomics/2.5 full_gwas_runs/step2_svat_lifestyles/{disease}/{file}",
                                     sep='\s+',
                                     usecols=['MarkerID', 'Allele2', 'Allele1', 'BETA', 'N_case', 'N_ctrl',
                                              'AF_Allele2', 'p.value', 'imputationInfo'])])

""" Combine """
data = pd.concat([geno, imputed])
data['MAF'] = 0.5 - abs(geno['AF_Allele2'] - 0.5)
data = (data
        .assign(N=data['N_case'] + data['N_ctrl'])
        # A1 = effect
        .rename(columns={'Allele2': 'A1', 'Allele1': 'A2', 'MarkerID': 'SNP', 'p.value': 'P', 'imputationInfo': 'INFO'})
        .drop_duplicates(subset=['SNP', 'A1', 'A2']))[['SNP', 'A1', 'A2', 'BETA', 'P', 'N', 'INFO']].copy()
# LDSC only works with strict single-nucleotide polymorphisms - i.e. remove all indels
data = data[(data['A1'].isin(['A', 'C', 'G', 'T'])) & (data['A2'].isin(['A', 'C', 'G', 'T']))].copy()

""" Export """
data.to_csv(f"~/ch2_genomics/2.10 whole_genome_genetic_correlations/ldsc/gwas_summary_stats_for_ldsc/collated/{disease}.tsv", sep='\t', index=False)
