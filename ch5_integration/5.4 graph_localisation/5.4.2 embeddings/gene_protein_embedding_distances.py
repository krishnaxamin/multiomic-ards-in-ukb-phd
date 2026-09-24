"""
Generate node2vec embeddings all genes in a tissue graph and calculate cosine similarities between each gene.
Then, determine whether GWAS genes and proteomics genes show greater within-set similarity than between-set similarity.
If so, this supports the hypothesis GWAS genes and proteomics genes are separated in the cellular network.
"""
import pandas as pd
import os
import sys
import networkx as nx
import scipy
import logging
import numpy as np
import torch
import torch_geometric
import argparse
import pickle
import math
import itertools

from pyprind import ProgBar
from typing import Dict, Union, List, Tuple, Set
from sklearn.metrics.pairwise import cosine_similarity

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

# digest Apocrita job name to give a log file
parser = argparse.ArgumentParser()
parser.add_argument('--task_index', type=int, default=1)
args = parser.parse_args()  # parse arguments

task_index = args.task_index
tissue_index = math.floor((task_index - 1) / 6)

proteomics_algos = ['firth', 'cox']
gene_edges_to_use_options = ['interacts_with', 'coexpresses_with', 'both']
options = list(itertools.product(proteomics_algos, gene_edges_to_use_options))
option = options[(task_index - 1) % 6]

ALGO = option[0]
gene_edges_to_use = option[1]
if gene_edges_to_use == 'both':
    gene_edges_used = ['interacts_with', 'coexpresses_with']
else:
    gene_edges_used = [gene_edges_to_use]

# select tissue
tissues = sorted([x.name for x in os.scandir(f"~/data/internal/knowledge_graph/tissue_graphs") if x.is_dir()])
TISSUE = tissues[tissue_index]

GENE_EDGES_USED_NAME = '_'.join([x.replace('_', '-') for x in gene_edges_used])
JOB_NAME = f"{GENE_EDGES_USED_NAME}_{ALGO}"
LOG_PATH = ''
LOGGING_PATH = (f"{LOG_PATH}/gene_protein_embedding_distances_{TISSUE}_{JOB_NAME}.log")

# In a module, __name__ = module's name in the Python package namespace
# __name__ is in the 'name' argument. logger with name foo is the parent of a logger with name foo.bar etc.
logger = logging.getLogger(f"GWAS-proteomics embeddings for {TISSUE}: {JOB_NAME}")
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
# DEFINE GENE AND PROTEIN SETS
#
########################################################################################################################


def get_proteomics_genes(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> List[str]:

    logger.info(f"Pulling genes associated to {disease_code} through proteomics ({ALGO}).")

    genes = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # there could be proteomic_firth, proteomics_cox, or magma edges between a gene and a disease
            if len(edges_data) > 3:
                logger.error(f"There are {len(edges_data)} edges between node {neighbour} and disease {disease_code}.")
                raise ValueError(
                    f"There are {len(edges_data)} edges between node {neighbour} and disease {disease_code}.")
            # iterate through the edges
            for edge_data in edges_data.values():
                # select edges which represent a proteomics Firth/Cox link
                if edge_data.get('edge_type') == f"proteomics_{ALGO}_associated_with":
                    genes.add(neighbour)

    logger.info(f"Pulling proteomics {ALGO} genes complete. Pulled {len(genes)} genes.")

    return list(genes)


def get_gwas_genes(full_graph: nx.MultiDiGraph, gene_gwas_relationship: List[str], disease_code: str = 'E11') -> List[str]:

    logger.info(f"Pulling genes associated to {disease_code} through GWAS.")

    genes = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # there could be proteomic_firth, proteomics_cox, or magma edges between a gene and a disease
            if len(edges_data) > 3:
                logger.error(f"There are {len(edges_data)} edges between node {neighbour} and disease {disease_code}.")
                raise ValueError(
                    f"There are {len(edges_data)} edges between node {neighbour} and disease {disease_code}.")
            # iterate through the edges (there may be multiple edges/edge types)
            for edge_data in edges_data.values():
                if 'variant-linked' in gene_gwas_relationship:
                    # select edges which represent a GWAS link
                    if edge_data.get('edge_type') == 'gwas_associated_with':
                        variant = neighbour
                        # now check genes linked to this variant - there can be many genes linked to the variant
                        for gene in full_graph.neighbors(variant):
                            variant_gene_edges = full_graph.get_edge_data(variant, gene)
                            if variant_gene_edges:
                                for variant_gene_edge in variant_gene_edges.values():
                                    if variant_gene_edge.get('edge_type') == 'proximal_to' and variant_gene_edge.get('feature_dist_to_proximal_gene') <= 2000:
                                        genes.add(gene)
                                # if e2[0].get('edge_type') == 'proximal_to':
                                #     genes.add(gene)
                if 'magma-linked' in gene_gwas_relationship:
                    # select edges which represent a MAGMA link
                    if edge_data.get('edge_type') == 'magma_associated_with':
                        genes.add(neighbour)

    logger.info(f"Pulling GWAS genes complete. Pulled {len(genes)} genes.")

    return list(genes)


########################################################################################################################
#
# GRAPH EMBEDDINGS: embedding, distance calcs, comparisons
#
########################################################################################################################


def create_homo_gene_subgraph(full_graph: nx.MultiDiGraph, inter_gene_edges: List[str]) -> nx.Graph:
    logger.info(f"Creating a homogeneous gene subgraph with edge(s): {' and '.join(inter_gene_edges)}.")

    # converting to sorted list to ensure ordering of the nodes in the subgraph is consistent between runs
    genes = sorted(list({n for n, d in full_graph.nodes(data=True) if d.get("type") == 'gene'}))
    logger.info(f"There are {len(genes)} genes in the full graph.")
    selected_edges = sorted([(u, v) for u, v, k, d in full_graph.edges(keys=True, data=True)
                             if u in genes and v in genes and d.get("edge_type") in inter_gene_edges])
    subgraph = nx.Graph()
    subgraph.add_nodes_from(genes)
    subgraph.add_edges_from(selected_edges)

    logger.info(f"Homogeneous gene subgraph created, with {len(subgraph.nodes())} nodes/genes. "
                f"Converting node labels to integer IDs.")

    # ordering argument forces a consistent mapping between the original node_ids and the integer IDs
    subgraph_int_ids = nx.relabel.convert_node_labels_to_integers(subgraph, label_attribute='original_id',
                                                                  ordering='sorted')

    logger.info('Conversion complete.')

    return subgraph_int_ids


def prep_networkx_for_pyg(
        networkx_graph: nx.MultiDiGraph, inter_gene_edges: List[str]
) -> Tuple[torch_geometric.data.Data, Dict[int, str]]:
    # subset homogeneous gene graph from full networkx graph, and re-ID so node IDs are consecutive from 0
    simple_graph = create_homo_gene_subgraph(full_graph=networkx_graph, inter_gene_edges=inter_gene_edges)
    logger.info(f"There are {len(simple_graph.nodes())} genes in the homogeneous subgraph.")

    pyg_to_nx_node_mapping = {new: data['original_id'] for new, data in simple_graph.nodes(data=True)}
    logger.info(f"pyg_to_nx_node_mapping covers {len(pyg_to_nx_node_mapping)} genes.")

    # convert networkx graph to pyg graph
    pyg_graph = torch_geometric.utils.from_networkx(simple_graph)

    return pyg_graph, pyg_to_nx_node_mapping


def node2vec_trainer(model: torch_geometric.nn.Node2Vec, loader: torch_geometric.loader.DataLoader, device: str,
                     optimiser: torch.optim.SparseAdam) -> float:
    model.train()
    total_loss = 0
    for pos_rw, neg_rw in loader:
        optimiser.zero_grad()
        loss = model.loss(pos_rw.to(device), neg_rw.to(device))
        loss.backward()
        optimiser.step()
        total_loss += loss.item()
    return total_loss / len(loader)


def train_node2vec(
        networkx_graph: nx.MultiDiGraph, inter_gene_edges: List[str]
) -> Dict[str, Union[np.ndarray, Dict[int, str], Dict[str, np.ndarray]]]:
    # subset homogeneous gene graph from full networkx graph, and re-ID so node IDs are consecutive from 0
    # convert networkx graph to pyg graph
    pyg_graph, pyg_to_nx_node_mapping = prep_networkx_for_pyg(networkx_graph=networkx_graph,
                                                              inter_gene_edges=inter_gene_edges)

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # parameters can be tuned - these are out-the-box parameters
    embedding_dim = 128
    walk_length = 20
    context_size = 10
    walks_per_node = 10

    # set up model
    model = torch_geometric.nn.Node2Vec(
        pyg_graph.edge_index,
        embedding_dim=embedding_dim,
        walk_length=walk_length,
        context_size=context_size,
        walks_per_node=walks_per_node,
        num_negative_samples=1,
        num_nodes=pyg_graph.num_nodes,
        p=1, q=1,
        sparse=True
    ).to(device)

    # set up loader and optimiser
    loader = model.loader(batch_size=128, shuffle=True, num_workers=0)
    optimiser = torch.optim.SparseAdam(list(model.parameters()), lr=0.01)

    # train over 50 epochs
    for epoch in range(1, 51):
        loss = node2vec_trainer(model=model, loader=loader, device=device, optimiser=optimiser)
        if epoch % 5 == 0:
            logger.info(f"Epoch {epoch}, Loss: {loss:.4f}")

    # get node embeddings
    embeddings = model.embedding.weight.detach().cpu().numpy()
    logger.info(f"Generated embeddings for {len(embeddings)} genes (using return from model).")
    # {nx_id: np.ndarray of shape (embedding_dim,) with embeddings}
    # nx_embeddings_dict = {pyg_graph.original_id[new]: emb for new, emb in enumerate(embeddings)}
    nx_embeddings_dict = {pyg_to_nx_node_mapping[new]: emb for new, emb in enumerate(embeddings)}

    logger.info(f"Generated embeddings for {len(nx_embeddings_dict)} genes (using return dict).")

    return {'embedding_nparray': embeddings, 'pyg_to_nx_map': pyg_to_nx_node_mapping, 'nx_embeddings': nx_embeddings_dict}


def calculate_embedding_similarities(
        embeddings_dict: Dict[str, np.ndarray], gwas_genes: List[str], proteomics_genes: List[str]) -> pd.DataFrame:
    gwas_proteomics_genes = list(set(gwas_genes) | set(proteomics_genes))
    relevant_embeddings_dict = {k: v for k, v in embeddings_dict.items() if k in gwas_proteomics_genes}

    ids = list(relevant_embeddings_dict.keys())

    genes_not_in_embeddings = list(set(gwas_proteomics_genes) - set(ids))
    if len(genes_not_in_embeddings) > 0:
        logger.error(f"{len(genes_not_in_embeddings)} GWAS-proteomics genes are not in the embeddings: "
                     f"{', '.join(genes_not_in_embeddings)}.")
        raise ValueError(f"{len(genes_not_in_embeddings)} GWAS-proteomics genes are not in the embeddings: "
                         f"{', '.join(genes_not_in_embeddings)}.")

    embeddings_mat = np.vstack([relevant_embeddings_dict[i] for i in ids])  # shape (n_nodes, dim)

    # cosine similarity matrix (n_nodes, n_nodes)
    cos_sim = cosine_similarity(embeddings_mat, embeddings_mat)

    cos_sim_df = pd.DataFrame(cos_sim, index=ids, columns=ids)

    return cos_sim_df


def calculate_rank_biserial_correlation(mann_whitney_u_statistic: float, n_group1: int, n_group2: int) -> float:
    r = ((2 * mann_whitney_u_statistic) / (n_group1 * n_group2)) - 1

    return r


def compare_node_group_similarities(similarity_df: pd.DataFrame, gwas: List[str], prot: List[str]) -> Dict[
    str, Dict[str, float]]:
    # used to check if: protein-protein and gwas-gwas similarities sig larger than protein-gwas distances

    np.fill_diagonal(similarity_df.values, np.nan)

    gwas_gwas_sims = similarity_df.loc[gwas, gwas].values.flatten()
    gwas_gwas_sims = gwas_gwas_sims[~np.isnan(gwas_gwas_sims)]
    prot_prot_sims = similarity_df.loc[prot, prot].values.flatten()
    prot_prot_sims = prot_prot_sims[~np.isnan(prot_prot_sims)]
    gwas_prot_sims = similarity_df.loc[gwas, prot].values.flatten()
    gwas_prot_sims = gwas_prot_sims[~np.isnan(gwas_prot_sims)]

    intra_gwas_vs_inter = scipy.stats.mannwhitneyu(gwas_gwas_sims, gwas_prot_sims, alternative='greater')
    intra_prot_vs_inter = scipy.stats.mannwhitneyu(prot_prot_sims, gwas_prot_sims, alternative='greater')
    intra_gwas_vs_intra_prot = scipy.stats.mannwhitneyu(gwas_gwas_sims, prot_prot_sims)

    return {'intra_gwas_vs_inter': {'U': intra_gwas_vs_inter.statistic, 'pval': intra_gwas_vs_inter.pvalue,
                                    'rank_biserial': calculate_rank_biserial_correlation(
                                        mann_whitney_u_statistic=intra_gwas_vs_inter.statistic,
                                        n_group1=len(gwas_gwas_sims), n_group2=len(gwas_prot_sims)
                                    )},
            'intra_prot_vs_inter': {'U': intra_prot_vs_inter.statistic, 'pval': intra_prot_vs_inter.pvalue,
                                    'rank_biserial': calculate_rank_biserial_correlation(
                                        mann_whitney_u_statistic=intra_prot_vs_inter.statistic,
                                        n_group1=len(prot_prot_sims), n_group2=len(gwas_prot_sims)
                                    )},
            'intra_gwas_vs_intra_prot': {'U': intra_gwas_vs_intra_prot.statistic,
                                         'pval': intra_gwas_vs_intra_prot.pvalue,
                                         'rank_biserial': calculate_rank_biserial_correlation(
                                             mann_whitney_u_statistic=intra_gwas_vs_intra_prot.statistic,
                                             n_group1=len(gwas_gwas_sims), n_group2=len(prot_prot_sims)
                                         )}}


def prep_embeddings_for_export(embeddings_dict: Dict[str, np.ndarray]) -> pd.DataFrame:
    first_emb = next(iter(embeddings_dict.values()))
    emb_dim = len(first_emb)

    df = pd.DataFrame.from_dict(
        embeddings_dict, orient='index', dtype=np.float32
    )
    df.columns = [f"dim_{i + 1}" for i in range(emb_dim)]
    df.index.name = 'gene'

    return df


########################################################################################################################
#
# MAIN EXECUTIONS
#
########################################################################################################################


def execution():

    tissue_graph = load_graph(tissue=TISSUE)

    # load pre-existing embeddings if present
    if os.path.exists(f"~/data/internal/knowledge_graph/tissue_graphs/{TISSUE}/n2v_data_{GENE_EDGES_USED_NAME}.pkl"):
        logger.info(f"Embeddings for {TISSUE} and {GENE_EDGES_USED_NAME} already exist. Pre-loading.")
        with open(f"~/data/internal/knowledge_graph/tissue_graphs/{TISSUE}/n2v_data_{GENE_EDGES_USED_NAME}.pkl", 'rb') as pickle_file:
            node2vec_training = pickle.load(pickle_file)
    else:
        # generate embeddings - embeddings valid across all diseases for a given tissue and gene-gene link set
        node2vec_training = train_node2vec(networkx_graph=tissue_graph,
                                           inter_gene_edges=['interacts_with', 'coexpresses_with'])

        # export node2vec_training
        logger.info(f"Exporting embeddings and PyG-to-NX mapping to '~/data/internal/knowledge_graph/tissue_graphs/{TISSUE}'.")
        with open(f"~/data/internal/knowledge_graph/tissue_graphs/{TISSUE}/n2v_data_{GENE_EDGES_USED_NAME}.pkl", 'wb') as pickle_file:
            pickle.dump(node2vec_training, pickle_file)

        # export embeddings
        embeddings_for_export = prep_embeddings_for_export(node2vec_training['nx_embeddings'])
        embeddings_for_export.to_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{TISSUE}/n2v_embeddings_{GENE_EDGES_USED_NAME}.csv")
        logger.info(f"Exporting complete.")

    # loop through all diseases and ways genes can be associated with disease via GWAS
    gene_gwas_relationships_list = [['variant-linked'], ['magma-linked'], ['variant-linked', 'magma-linked']]
    single_tissue_results_list = []
    disease_info = pd.read_csv(f"~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv")
    bar = ProgBar(len(disease_info), stream=sys.stdout, title=f"Analysing all diseases for {TISSUE}")
    for disease in sorted(disease_info.icd10_three_letter.to_list()):

        logger.info(f"-----RUNNING FOR {disease}-----")
        for gene_gwas_relationship in gene_gwas_relationships_list:
            logger.info(f"--------RUNNING FOR {' & '.join(gene_gwas_relationship)} GWAS gene-disease relationships.")

            job_name_with_gene_gwas = f"{JOB_NAME}_{'_'.join(gene_gwas_relationship)}"

            proteomics_genes = get_proteomics_genes(full_graph=tissue_graph, disease_code=disease)
            gwas_genes = get_gwas_genes(full_graph=tissue_graph, gene_gwas_relationship=gene_gwas_relationship,
                                        disease_code=disease)

            # there must be both proteomics-associated and GWAS-associated genes for the analysis to be worth running
            # detect gene counts for both gene sets before performing the skip
            if len(proteomics_genes) == 0:
                logger.info(f"There are no proteomics ({ALGO}) genes for {disease}. Skipping to next disease.")
            if len(gwas_genes) == 0:
                logger.info(f"There are no {' or '.join(gene_gwas_relationship)} GWAS genes for {disease}. Skipping to next disease.")
            if len(proteomics_genes) * len(gwas_genes) == 0:
                continue

            # test embedding separation between gwas-/proteomic-set genes
            embedding_similarities = calculate_embedding_similarities(embeddings_dict=node2vec_training['nx_embeddings'],
                                                                      gwas_genes=gwas_genes, proteomics_genes=proteomics_genes)
            gwas_proteomic_embedding_sep = compare_node_group_similarities(similarity_df=embedding_similarities,
                                                                           gwas=gwas_genes, prot=proteomics_genes)
            logger.info(f"Embedding separation test complete: real-world data.")

            logger.info(f"Comparing similarities between node2vec embeddings of GWAS- and proteomics-genes: \n"
                        f"- intra-GWAS vs inter-GWAS-proteomics 'greater' pval: {gwas_proteomic_embedding_sep['intra_gwas_vs_inter']['pval']} \n"
                        f"- intra-GWAS vs inter-GWAS-proteomics rank-biserial: {gwas_proteomic_embedding_sep['intra_gwas_vs_inter']['rank_biserial']} \n"
                        f"- intra-proteomics vs inter-GWAS-proteomics 'greater' pval: {gwas_proteomic_embedding_sep['intra_prot_vs_inter']['pval']} \n"
                        f"- intra-proteomics vs inter-GWAS-proteomics rank-biserial: {gwas_proteomic_embedding_sep['intra_prot_vs_inter']['rank_biserial']} \n"
                        f"- intra-GWAS vs intra-proteomics 'two-sided' pval: {gwas_proteomic_embedding_sep['intra_gwas_vs_intra_prot']['pval']} \n"
                        f"- intra-proteomics vs inter-GWAS-proteomics rank-biserial: {gwas_proteomic_embedding_sep['intra_gwas_vs_intra_prot']['rank_biserial']}"
                        f"")

            single_tissue_results_list.append({'disease': disease, 'tissue': TISSUE,
                                               'gwas_gene_disease_relationship': '_'.join(gene_gwas_relationship),
                                               'intra_gwas_vs_inter_pval':
                                                   gwas_proteomic_embedding_sep['intra_gwas_vs_inter']['pval'],
                                               'intra_gwas_vs_inter_rank_biserial':
                                                   gwas_proteomic_embedding_sep['intra_gwas_vs_inter']['rank_biserial'],
                                               'intra_proteomics_vs_inter_pval':
                                                   gwas_proteomic_embedding_sep['intra_prot_vs_inter']['pval'],
                                               'intra_proteomics_vs_inter_rank_biserial':
                                                   gwas_proteomic_embedding_sep['intra_prot_vs_inter']['rank_biserial'],
                                               'intra_gwas_vs_intra_prot_pval':
                                                   gwas_proteomic_embedding_sep['intra_gwas_vs_intra_prot']['pval'],
                                               'intra_gwas_vs_intra_prot_rank_biserial':
                                                   gwas_proteomic_embedding_sep['intra_gwas_vs_intra_prot']['rank_biserial'],
                                               })

            bar.update()

    single_tissue_results_df = pd.DataFrame(single_tissue_results_list)
    single_tissue_results_df.to_csv(f"~/data/internal/knowledge_graph/gene_protein_distances/embeddings/{TISSUE}/embedding_comparisons_{JOB_NAME}.csv",
                                    index=False)

    logger.info(f"Execution for {TISSUE} complete.")


if __name__ == '__main__':
    execution()
