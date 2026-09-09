""" Extract covariate information from the full genomics cohort. """
import pandas as pd

data = pd.read_csv('~/data/external/genomics/array-genotyping_phenotype_file.csv')

covariate_fields = ['eid', 'p21003_i0', 'p31', 'p52', 'p34', 'p53_i0', 'p54_i0', 'p22000', 'p22007', 'p22008']

genomics_covariates = data[covariate_fields]

genomics_covariates.to_csv('~/data/external/internal/genomics_covars.csv', index=False)
