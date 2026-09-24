# Code
- `gene_protein_embedding_distances.py`: Generate node2vec embeddings all genes in a tissue graph and calculate cosine similarities between each gene. Then, determine whether GWAS genes and proteomics genes show greater within-set similarity than between-set similarity.
- `gene_protein_embedding_distances.sh`: Shell script to execute `gene_protein_embedding_distances.py`.
- `plot_embedding_distance_results.py`: Plot results of `gene_protein_embedding_distances.py`.
- `plot_embedding_distance_biserials.py`: Plot proportions of tissue-disease instances that have positive or negative rank biserials, resulting from `gene_protein_embedding_distances.py`.
- `olink_cossims_vs_all.py`: Using gene embeddings to test whether Olink proteins are particularly concentrated within cellular networks.
- `olink_cossims_vs_all.sh`: Shell script to execute `olink_cossims_vs_all.py`. 