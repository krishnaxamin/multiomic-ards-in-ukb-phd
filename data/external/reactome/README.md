# Data present
- `reactome_id_ancestor_relations.csv`: Get all ancestors for each Reactome annotation, using `~/utils/reactome_tree_extraction.py`. Required because `ReactomePathwaysRelation.txt` lists only the immediate parent-chil relationships.

# Data absent
- `NCBI2Reactome_PE_All_Levels.txt`: Maps genes identified by NCBI IDs to their Reactome annotations at all levels of the pathway hierarchy, not just the lowest level. Downloaded from https://download.reactome.org/92/NCBI2Reactome_PE_All_Levels.txt on 23 March, 2025.
- `ReactomePathwaysRelation.txt`: Reactome pathways hierarchy relationships. Downloaded from https://download.reactome.org/92/ReactomePathwaysRelation.txt on 25 March, 2025.
- `UniProt2Reactome_PE_All_Levels.txt`: Maps genes identified by UniProt IDs to their Reactome annotations at all levels of the pathway hierarchy, not just the lowest level. Downloaded from https://download.reactome.org/92/UniProt2Reactome_PE_All_Levels.txt on 25 March, 2025.
- `Ensembl2Reactome_PE_All_Levels.txt`: Maps genes identified by Ensembl IDs (ENSG IDs) to their Reactome annotations at all levels of the pathway hierarchy, not just the lowest level. Downloaded from https://download.reactome.org/92/Ensembl2Reactome_PE_All_Levels.txt on 25 March, 2025.