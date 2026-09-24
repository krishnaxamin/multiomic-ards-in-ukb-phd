""" Construct file mapping gene symbols to UniProt codes to UKB fields for the Olink proteins. """
from pandas import read_csv, DataFrame, concat


def expand_row(row, col1, col2, col3):
    c1, c2, c3 = row[col1], row[col2], row[col3]

    s1, s2 = c1.split("_"), c2.split("_")

    if "_" in c1 and "_" in c2 and len(s1) == len(s2):
        # Case 1: same number of underscores → pair positionally
        return DataFrame({col1: s1, col2: s2, col3: [c3]*len(s1)})

    elif "_" not in c1 and "_" not in c2:
        # Case 2: no underscores → keep as is
        return DataFrame([row])

    else:
        # Case 3: only one column has underscores → expand only that one
        if "_" in c1:
            return DataFrame({col1: s1, col2: [c2]*len(s1), col3: [c3]*len(s1)})
        else:
            return DataFrame({col1: [c1]*len(s2), col2: s2, col3: [c3]*len(s2)})



uniprot_genesymbol = read_csv('~/data/external/proteomics/olink_assay.dat', sep='\t')[['Assay', 'UniProt']].rename(columns={'Assay': 'gene_symbol', 'UniProt': 'uniprot'})
uniprot_genesymbol['protein_field'] = uniprot_genesymbol.gene_symbol.str.lower()
uniprot_genesymbol['protein_field'] = uniprot_genesymbol.protein_field.replace({'hla-dra': 'hla_dra', 'hla-a': 'hla_a', 'hla-e': 'hla_e', 'ervv-1': 'ervv_1'})
uniprot_genesymbol_expanded = concat([expand_row(row, col1='gene_symbol', col2='uniprot', col3='protein_field')
                                      for _, row in uniprot_genesymbol.iterrows()], ignore_index=True)
uniprot_genesymbol_expanded.to_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv', index=False)
