# Code
- `calculate_semantic_similarity.py`: Calculate pathway-based semantic similarity between ARDs.
- `clustering.py`: Cluster ARDs using Louvain clustering on a graph connecting ARDs using semantic similarity.
- `identifying_driver_pathways.py`: Identify pathways driving semantic similarity between ARDs, and identify sets of ARDs that share the same 'complete connector' driver pathways.
- `by_magma_sensitivity_analyses.py`: Sensitivity analyses to explore effect of using different (1) top-N pathways in calculating semantic similarity (2) top-N MICAs in identifying driver pathways.