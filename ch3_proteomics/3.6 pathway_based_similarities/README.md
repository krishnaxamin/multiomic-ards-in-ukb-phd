# Code
- `calculate_semantic_similarity.py`: Calculate pathway-based semantic similarity between ARDs.
- `clustering.py`: Cluster ARDs using Louvain clustering on a graph connecting ARDs using semantic similarity.
- `identifying_driver_pathways.py`: Identify pathways driving semantic similarity between ARDs, and identify sets of ARDs that share the same 'complete connector' driver pathways.
- `firth_vs_cox_semsim.py`: Calculate semantic similarity between ARDs and results from different regressions, i.e., disease A from Firth x disease B from Cox. Find driver pathways. Investigate prevalent-incident relationships involving incident ICD-10 chapter 13 ARDs.