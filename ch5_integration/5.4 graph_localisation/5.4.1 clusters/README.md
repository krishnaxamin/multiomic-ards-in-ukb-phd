# Code
- `gene_protein_in_diff_clusters.py`: Detect clusters in the graph, and test whether GWAS genes and proteomics genes are enriched in each of the clusters. Then, test whether there is more co-occurrence of GWAS and proteomic enrichment in the same clusters than expected, for each tissue-disease combination. Overall, tests the hypothesis that GWAS gene and proteomics genes are separated in the cellular network.
- `gene_protein_in_diff_clusters.sh`: Shell script to execute `gene_protein_in_diff_clusters.py`.
- `plot_diff_cluster_results.py`: Plot the enrichment distribution of clusters in each ARD-tissue combination. Enrichment distribution = proportion of clusters enriched in either GWAS or proteomics, or both.
- `generate_cluster_enrichment_stats.py`: Generate counts of clusters enriched in GWAS-only, proteomics-only, both, or neither.
- `plot_n_clusters_per_tissue.py`: Plot the number of clusters for each tissue/gene-gene link combination. 