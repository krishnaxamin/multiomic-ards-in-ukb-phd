"""
Script to make pheno files for GWAS. Includes: FID/IID (or both), disease phenotype, covariates.
"""

import sys
from pandas import read_csv, merge, DataFrame
from pyprind import ProgBar

non_lifestyle_covariates = read_csv('~/data/internal/genomics/genomics_covars.csv')
disease_fields = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'].to_list()
all_phenotypes = read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv')

non_lifestyle_covariates = non_lifestyle_covariates.replace('Stockport (pilot)', 'Stockport')
covariates = non_lifestyle_covariates[['eid', 'p21003_i0', 'p31', 'p54_i0', 'p22000', 'p22008']].copy()

# categorical covariates for GCTA
gcta_fixed_covars = non_lifestyle_covariates.rename(columns={'eid': 'IID'})
gcta_fixed_covars['FID'] = gcta_fixed_covars['IID']  # for a person, FID = IID
gcta_fixed_covars = gcta_fixed_covars[[gcta_fixed_covars.columns.to_list()[-1]] + gcta_fixed_covars.columns.to_list()[:-1]]
gcta_fixed_covars = gcta_fixed_covars.drop(['p21003_i0', 'p22007', 'p22008', 'p22000'], axis=1)  # p22000 is categorical, but has too many categories so has to be considred as quantitative
gcta_fixed_covars.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/categorical_covars_GCTA.tsv', sep='\t', index=False)

# REGENIE, SAIGE, GCTA quantitative covariates and phenotype
bar = ProgBar(len(disease_fields), stream=sys.stdout, title='Making pheno files')
for disease in disease_fields:
    # i = 30
    # disease = 'p131036'
    pca_results = read_csv('~/data/internal/genomics/pca/' + disease + '/' + disease + '_pca.eigenvec', delim_whitespace=True).drop('#FID', axis=1)
    phenotypes = all_phenotypes[['eid', disease]]
    covariates_pca = merge(covariates, pca_results, left_on='eid', right_on='IID').drop('IID', axis=1)

    # GCTA quantitative covars
    gcta_quant_covars = covariates_pca.drop(['p31', 'p54_i0', 'p22007', 'p22008'], axis=1)
    gcta_quant_covars = gcta_quant_covars.rename(columns={'eid': 'IID'})
    gcta_quant_covars['FID'] = gcta_quant_covars['IID']  # for a person, FID = IID
    gcta_quant_covars = gcta_quant_covars[
        [gcta_quant_covars.columns.to_list()[-1]] + gcta_quant_covars.columns.to_list()[:-1]]
    batch_vars_to_integer_map = DataFrame({'p22000': list(set(gcta_quant_covars['p22000'])),
                                           'batch_integer': range(1, len(set(gcta_quant_covars['p22000'])) + 1)})
    gcta_quant_covars = merge(gcta_quant_covars, batch_vars_to_integer_map, on='p22000')
    gcta_quant_covars.drop('p22000', axis=1, inplace=True)
    gcta_quant_covars.rename(columns={'batch_integer': 'p22000'}, inplace=True)
    gcta_quant_covars.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/' + disease + '_quant_covars_GCTA.tsv', sep='\t', index=False)

    # GCTA phenotypes (also valid for REGENIE)
    gcta_pheno = phenotypes.rename(columns={'eid': 'IID'})
    gcta_pheno['FID'] = gcta_pheno['IID']  # for a person, FID = IID
    gcta_pheno = gcta_pheno[[gcta_pheno.columns.to_list()[-1]] + gcta_pheno.columns.to_list()[:-1]]
    gcta_pheno.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/fastGWA-GLMM/' + disease + '_pheno-file_GCTA.tsv', sep='\t', index=False)
    gcta_pheno.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/' + disease + '_pheno-file_REGENIE.tsv', sep='\t', index=False)

    # REGENIE, SAIGE pheno files (phenotypes + all covars) (BOLT file w/o phenotype = REGENIE covar file)
    pheno_file = merge(covariates_pca, phenotypes, on='eid')
    pheno_file.rename(columns={'eid': 'IID'}, inplace=True)
    pheno_file.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/SAIGE/' + disease + '_pheno-file_SAIGE.tsv', sep='\t', index=False)
    pheno_file['FID'] = pheno_file['IID']  # for a person, FID = IID
    pheno_file = pheno_file[[pheno_file.columns.to_list()[-1]] + pheno_file.columns.to_list()[:-1]]
    pheno_file = pheno_file.drop(['p22007', 'p22008'], axis=1)
    pheno_file.drop(disease, axis=1, inplace=True)
    pheno_file.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/REGENIE/' + disease + '_covars_REGENIE.tsv', sep='\t', index=False)
    bar.update()


""" Pheno files, including lifestyle factors """
# to be run on Apocrita
technical_covariates = read_csv('~/data/internal/genomics/genomics_covars.csv')
lifestyle_covariates = read_csv('~/data/internal/genomics/genomics_lifestyles.csv')
disease_fields = read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'].to_list()
all_phenotypes = read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv')

lifestyle_covariates = lifestyle_covariates.replace('Stockport (pilot)', 'Stockport')
covariates = technical_covariates.merge(lifestyle_covariates, on='eid')
covariates = covariates[['eid', 'p21003_i0', 'p31', 'p54_i0', 'p23099_i0', 'p189', 'p884_i0', 'p1160_i0',
                         'p1289_i0', 'p1299_i0', 'p1309_i0', 'p1319_i0', 'p1438_i0', 'p1458_i0', 'p1528_i0']].copy()

# BOLT-LMM, SAIGE, GCTA quantitative covariates and phenotype
bar = ProgBar(len(disease_fields), stream=sys.stdout, title='Making pheno files')
for disease in disease_fields:
    # i = 30
    # disease = 'p131036'
    pca_results = read_csv('~/data/internal/genomics/pca/' + disease + '/' + disease + '_pca.eigenvec', delim_whitespace=True).drop('#FID', axis=1)
    phenotypes = all_phenotypes[['eid', disease]]
    # non_lifestyle_covariates_pca = merge(non_lifestyle_covariates, pca_results, left_on='eid', right_on='IID').drop('IID', axis=1)
    covariates_pca = merge(covariates, pca_results, left_on='eid', right_on='IID').drop('IID', axis=1)

    # BOLT-LMM, SAIGE pheno files (phenotypes + all covars) (BOLT file w/o phenotype = REGENIE covar file)
    pheno_file = merge(covariates_pca, phenotypes, on='eid')
    pheno_file.rename(columns={'eid': 'IID'}, inplace=True)
    pheno_file.to_csv('~/data/internal/genomics/pheno_covar_files_for_gwas/SAIGE/' + disease + '_pheno-file_SAIGE_lifestyles.tsv', sep='\t', index=False)
    bar.update()
