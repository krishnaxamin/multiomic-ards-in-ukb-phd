"""
Using gene embeddings to test whether Olink proteins are concentrated within cellular networks.
Specifically, compares Olink-Olink cosine similarities to cosine similarities between all other background gene pairs.
 Background gene pairs do not include Olink proteins as a constituent.
"""
import pandas as pd
import os
import networkx as nx
import scipy
import logging
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import umap
import sklearn

from typing import Dict, List, Set
from sklearn.metrics.pairwise import cosine_similarity

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

LOG_PATH = ''
LOGGING_PATH = (f"{LOG_PATH}/olink_cossims_vs_all.log")

# In a module, __name__ = module's name in the Python package namespace
# __name__ is in the 'name' argument. logger with name foo is the parent of a logger with name foo.bar etc.
logger = logging.getLogger(f"Embedding cossims - intra-Olink vs all others")
logger.setLevel(logging.INFO)
logger.propagate = False

file_handler = logging.FileHandler(LOGGING_PATH, mode='w')
file_handler.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s - %(name)s. %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)


########################################################################################################################
#
# LOAD GRAPH
#
########################################################################################################################


def load_graph(tissue: str = 'Adipose_Subcutaneous') -> nx.MultiDiGraph:
    logger.info(f"Loading graph for {tissue}.")

    graph = nx.read_graphml(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/graph.graphml")

    logger.info(f"Graph for {tissue} loaded.")

    return graph


########################################################################################################################
#
# DEFINE OLINK SETS
#
########################################################################################################################


def get_olink_proteins() -> List[str]:
    logger.info(f"Pulling proteins on the Olink panel.")

    olink_proteins = pd.read_csv(f"~/data/internal/proteomics/olink_field_uniprot_ensg_genesymbol_mapping.csv")
    high_missingness_proteins = pd.read_csv(f"~/data/internal/proteomics/high_missingness_protein_fields.txt")

    olink_post_qc_proteins = olink_proteins[
        ~olink_proteins.protein_field.isin(high_missingness_proteins.protein_field)].copy()

    olink_ensembl = sorted(list(olink_post_qc_proteins.ensg_id.unique()))

    logger.info(f"{len(olink_proteins)} proteins on the Olink panel. {len(olink_ensembl)} pass our QC.")

    return olink_ensembl


def get_graph_olink_genes(full_graph: nx.MultiDiGraph, olink_genes: List[str]) -> List[str]:

    # get the Olink genes that are expressed in this tissue
    graph_olink_genes = [n for n in full_graph.nodes() if n in olink_genes]

    return graph_olink_genes


########################################################################################################################
#
# GRAPH EMBEDDINGS: embedding, distance calcs, comparisons
#
########################################################################################################################


def load_embeddings(tissue: str = 'Adipose_Subcutaneous'):
    logger.info(f"Loading in embeddings for {tissue} and all gene-gene links.")
    with open(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/n2v_data_interacts-with_coexpresses-with.pkl",
              'rb') as pickle_file:
        node2vec_training = pickle.load(pickle_file)

    return node2vec_training


def calculate_embedding_similarities(embeddings_dict: Dict[str, np.ndarray]) -> pd.DataFrame:
    logger.info(f"Calculating embedding cosine similarities.")

    ids = list(embeddings_dict.keys())

    embeddings_mat = np.vstack([embeddings_dict[i] for i in ids])  # shape (n_nodes, dim)

    # cosine similarity matrix (n_nodes, n_nodes)
    cos_sim = cosine_similarity(embeddings_mat, embeddings_mat)

    cos_sim_df = pd.DataFrame(cos_sim, index=ids, columns=ids)

    return cos_sim_df


def olink_olink_vs_olink_nonolink(olink_genes: List[str], cossims: pd.DataFrame):
    # this is the same logic as intra_proteomics_vs_inter in gene_protein_embedding_distances.py
    logger.info(f"Comparing cossims between Olink genes to the rest of background.")
    np.fill_diagonal(cossims.values, np.nan)

    olink_cossims = cossims.loc[olink_genes, olink_genes].values.flatten()
    olink_cossims = olink_cossims[~np.isnan(olink_cossims)]

    non_olink_cossims = np.concatenate((cossims.loc[olink_genes, :].drop(columns=olink_genes).values.flatten(),
                                        cossims.loc[:, olink_genes].drop(index=olink_genes).values.flatten()))
    non_olink_cossims = non_olink_cossims[~np.isnan(non_olink_cossims)]

    two_sided_pval = scipy.stats.mannwhitneyu(olink_cossims, non_olink_cossims).pvalue
    if two_sided_pval < 0.05:
        greater_pval = scipy.stats.mannwhitneyu(olink_cossims, non_olink_cossims, alternative='greater').pvalue
        less_pval = scipy.stats.mannwhitneyu(olink_cossims, non_olink_cossims, alternative='less').pvalue
    else:
        greater_pval = np.nan
        less_pval = np.nan

    return {'two_sided_pval': two_sided_pval, 'greater_pval': greater_pval, 'less_pval': less_pval}


def olink_nonolink_vs_nonolink_nonolink(olink_genes: List[str], cossims: pd.DataFrame):
    # this is the same logic as intra_proteomics_vs_inter in gene_protein_embedding_distances.py
    logger.info(f"Comparing cossims between Olink genes to the rest of background.")
    np.fill_diagonal(cossims.values, np.nan)

    olink_nonolink_cossims = np.concatenate((cossims.loc[olink_genes, :].drop(columns=olink_genes).values.flatten(),
                                             cossims.loc[:, olink_genes].drop(index=olink_genes).values.flatten()))
    olink_nonolink_cossims = olink_nonolink_cossims[~np.isnan(olink_nonolink_cossims)]

    nonolink_nonolink_cossims = cossims.drop(columns=olink_genes, index=olink_genes).values.flatten()
    nonolink_nonolink_cossims = nonolink_nonolink_cossims[~np.isnan(nonolink_nonolink_cossims)]

    two_sided_pval = scipy.stats.mannwhitneyu(olink_nonolink_cossims, nonolink_nonolink_cossims).pvalue
    if two_sided_pval < 0.05:
        greater_pval = scipy.stats.mannwhitneyu(olink_nonolink_cossims, nonolink_nonolink_cossims,
                                                alternative='greater').pvalue
        less_pval = scipy.stats.mannwhitneyu(olink_nonolink_cossims, nonolink_nonolink_cossims,
                                             alternative='less').pvalue
    else:
        greater_pval = np.nan
        less_pval = np.nan

    return {'two_sided_pval': two_sided_pval, 'greater_pval': greater_pval, 'less_pval': less_pval}


########################################################################################################################
#
# MAIN EXECUTIONS
#
########################################################################################################################


def single_tissue_execution(tissue: str, olink_genes: List[str]) -> Dict[str, Dict[str, float]]:
    tissue_graph = load_graph(tissue=tissue)

    olink_graph_genes = get_graph_olink_genes(tissue_graph, olink_genes=olink_genes)
    logger.info(f"{len(olink_graph_genes)} are expressed in {tissue}.")

    embeddings_from_pickle = load_embeddings(tissue)

    cossims = calculate_embedding_similarities(embeddings_dict=embeddings_from_pickle['nx_embeddings'])

    olink_olink_vs_olink_nonolink_pvals = olink_olink_vs_olink_nonolink(olink_genes=olink_graph_genes,
                                                                        cossims=cossims)
    olink_nonolink_vs_nonolink_nonolink_pvals = olink_nonolink_vs_nonolink_nonolink(
        olink_genes=olink_graph_genes,
        cossims=cossims)

    logger.info(f"Olink-Olink vs Olink-nonOlink: {olink_olink_vs_olink_nonolink_pvals}")
    logger.info(f"Olink-nonOlink vs nonOlink-nonOlink: {olink_nonolink_vs_nonolink_nonolink_pvals}")

    return {'olink_olink_vs_olink_nonolink_pvals': olink_olink_vs_olink_nonolink_pvals,
            'olink_nonolink_vs_nonolink_nonolink_pvals': olink_nonolink_vs_nonolink_nonolink_pvals}


def execution():
    tissues = sorted([x.name for x in os.scandir(f"~/data/internal/knowledge_graph/tissue_graphs") if x.is_dir()])

    # get all the genes/proteins on the Olink panel that pass our QC
    olink_genes = get_olink_proteins()

    olink_vs_all_all_tissues_list = []
    for tissue in tissues:
        logger.info(f"------Executing for {tissue}.------")
        single_tissue_results = single_tissue_execution(tissue=tissue, olink_genes=olink_genes)
        olink_vs_all_all_tissues_list.append({'tissue': tissue} |
                                             {f"olink_olink_vs_olink_nonolink_{k}": v for k, v in
                                              single_tissue_results['olink_olink_vs_olink_nonolink_pvals'].items()} |
                                             {f"olink_nonolink_vs_nonolink_nonolink_{k}": v for k, v in
                                              single_tissue_results[
                                                  'olink_nonolink_vs_nonolink_nonolink_pvals'].items()})

    olink_vs_all_all_tissues = pd.DataFrame(olink_vs_all_all_tissues_list)

    logger.info(f"Exporting results to '~/data/internal/knowledge_graph/gene_protein_distances/embeddings/"
                f"olink_vs_all_embedding_cossims.csv'")
    olink_vs_all_all_tissues.to_csv(f"~/data/internal/knowledge_graph/gene_protein_distances/embeddings/"
                                    f"olink_vs_all_embedding_cossims.csv", index=False)

    logger.info(f"Export and execution complete.")


if __name__ == '__main__':
    execution()
