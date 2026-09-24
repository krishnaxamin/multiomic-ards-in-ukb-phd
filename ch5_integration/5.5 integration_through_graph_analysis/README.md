# Code
- `determine_gwas_protoemics_disease_analysis_order.py`: Compare age-of-onsets for certain ARDs to determine relative 'order', i.e. whose genomics results are used, and whose proteomics results are used, given a set of ARDs. 
- `gene_protein_tissue_specific_links.py`: Integrate genomics and proteomics results using graph methods on the knowledge graph.
- `gene_protein_tissue_specific_links.sh`: Shell script to execute `gene_protein_tissue_specific_links.py`.
- `analyse_link_results.py`: Analyse results from `gene_protein_tissue_specific_links.py`. All results in R5.5 come from this script.