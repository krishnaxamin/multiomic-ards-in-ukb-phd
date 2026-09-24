""" Generate pairs of ARDs for use in HDL. """
from pandas import read_csv, DataFrame
from itertools import combinations

disease_fields = list(read_csv('~/data/interal/phenotype_coding/pan_ukbb_eur_qu10_50_disease_fields.txt')['disease_field'])
disease_pairs = [x[0] + '-' + x[1] for x in list(combinations(disease_fields, 2))]
DataFrame({'disease_pair': disease_pairs}).to_csv(
    '~/ch2_genomics/2.10 whole_genome_genetic_correlations/disease_pairs.txt', index=False)
