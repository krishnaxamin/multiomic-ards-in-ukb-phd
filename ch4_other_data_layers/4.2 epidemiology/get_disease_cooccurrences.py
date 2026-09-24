""" Get disease co-occurrence counts from phenotype files. """
from pandas import read_csv, DataFrame
from pyprind import ProgBar
import itertools
import sys


def get_cooccurrences(phenotypes):
    cooccurrences_matrix = DataFrame(index=list(phenotypes.columns)[1:],
                                     columns=list(phenotypes.columns)[1:])
    cooccurrences_matrix[cooccurrences_matrix.isna()] = 0

    bar = ProgBar(len(phenotypes), stream=sys.stdout)
    for i in range(len(phenotypes)):
        eid_phenotypes = phenotypes.iloc[i, 1:].copy()
        if eid_phenotypes.sum() > 1:
            indices_for_diseases = list(eid_phenotypes[eid_phenotypes == 1].index)
            permutations = list(itertools.permutations(indices_for_diseases, 2))
            for tup in permutations:
                # print(tup)
                cooccurrences_matrix.at[tup[0], tup[1]] = cooccurrences_matrix.at[tup[0], tup[1]] + 1
        bar.update()

    cooccurrences_matrix['disease1'] = list(cooccurrences_matrix.columns)
    cooccurrences_long = cooccurrences_matrix.melt(id_vars='disease1',
                                                   var_name='disease2',
                                                   value_name='cooccurrence_count')
    cooccurrences_long = cooccurrences_long[
        cooccurrences_long.disease1 != cooccurrences_long.disease2]

    return cooccurrences_long


""" Genomics """
genomics_phenotypes = read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_phenotypes.csv')

genomics_cooccurrences_matrix = DataFrame(index=list(genomics_phenotypes.columns)[1:],
                                          columns=list(genomics_phenotypes.columns)[1:])
genomics_cooccurrences_matrix[genomics_cooccurrences_matrix.isna()] = 0

bar = ProgBar(len(genomics_phenotypes), stream=sys.stdout)
for i in range(len(genomics_phenotypes)):
    eid_phenotypes = genomics_phenotypes.iloc[i, 1:].copy()
    if eid_phenotypes.sum() > 1:
        indices_for_diseases = list(eid_phenotypes[eid_phenotypes == 1].index)
        permutations = list(itertools.permutations(indices_for_diseases, 2))
        for tup in permutations:
            # print(tup)
            genomics_cooccurrences_matrix.at[tup[0], tup[1]] = genomics_cooccurrences_matrix.at[tup[0], tup[1]] + 1
    bar.update()

genomics_cooccurrences_matrix['disease1'] = list(genomics_cooccurrences_matrix.columns)
genomics_cooccurrences_long = genomics_cooccurrences_matrix.melt(id_vars='disease1',
                                                                 var_name='disease2',
                                                                 value_name='cooccurrence_count')
genomics_cooccurrences_long = genomics_cooccurrences_long[
    genomics_cooccurrences_long.disease1 != genomics_cooccurrences_long.disease2]
genomics_cooccurrences_long.to_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences.csv',
                                   index=False)

""" Proteomics, phenotypes prior to blood taking """
proteomics_prior_phenotypes = read_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_prior_phenotypes.csv')
proteomics_prior_cooccurrences = get_cooccurrences(proteomics_prior_phenotypes)
proteomics_prior_cooccurrences.to_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences.csv',
    index=False)

""" Proteomics, full cohort """
proteomics_phenotypes = read_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_phenotypes.csv')
proteomics_cooccurrences = get_cooccurrences(proteomics_phenotypes)
proteomics_cooccurrences.to_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences.csv',
                                index=False)
