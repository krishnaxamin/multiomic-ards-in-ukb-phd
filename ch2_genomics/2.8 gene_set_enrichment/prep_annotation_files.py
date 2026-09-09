""" Prep gene annotation files for use in MAGMA step 3. """
from pandas import read_csv, merge
from utils.enrichment_analyses import slim_annotations

# data downloaded: 10 June 2025
# https://download.reactome.org/92/NCBI2Reactome_PE_All_Levels.txt
reactome = read_csv('~/data/external/reactome/NCBI2Reactome_PE_All_Levels.txt', sep='\t',
                    header=None)
reactome.columns = ['entrez_id', 'reactome_id', 'reactome_name', 'reactome_pathway_id', 'reactome_pathway_link',
                    'reactome_pathway_name', 'evidence', 'species']
reactome_human_tas = reactome[(reactome['species'] == 'Homo sapiens') & (reactome['evidence'] == 'TAS')].copy()
reactome_human_tas_gene_annotation = reactome_human_tas[['entrez_id', 'reactome_pathway_id']].drop_duplicates(
    ignore_index=True).copy()
# slim to reduce too-general and too-specific annotations
reactome_human_tas_gene_annotation = slim_annotations(reactome_human_tas_gene_annotation,
                                                      annotation_col_id='reactome_pathway_id',
                                                      entity_col_id='entrez_id')
reactome_human_tas_gene_annotation.to_csv('~/data/internal/magma/reactome_annotations_for_magma.tsv', sep='\t',
                                          header=False, index=False)
