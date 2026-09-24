# Code
- `enrichment_analyses.py`: Functions to conduct, and relevant for, annotation enrichment analyses. 
- `generate_diease_pairs.py`: Generate **_C_**(68,2) = 2278 disease pairs for use in HDL, LDSC, and LCV.
- `get_genotyped_ref_alt.py`: Use UK Biobank Resource 1955 (https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955) to get REF/ALT alleles for the
genotyped variants in the UK Biobank.
- `get_number_of_variants_tested.py`: Get total number of variants tested in the GWAS, per-disease and per-disease pair.
- `mainali_alpha.py`: Means to implement Mainali et al 2022's alpha co-occurrence metric in Python.
- `olink_field_uniprot_genesymbol_mapping.py`: Construct file mapping gene symbols to UniProt codes to UKB fields for the Olink proteins. Outputs `~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv`.
- `reactome_tree_extraction`: Get all ancestors for all Reactome annotations. Outputs `~/data/external/reactome/reactome/reactome_id_ancestor_relations.csv`.
- `semantic_similarity.py`: Functions for pathway-based semantic similarity work.
- `significance_labelling.py`: Function to label a dataframe with nominal, FDR, FDR-permissive and Bonferroni significance levels.