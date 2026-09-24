import pandas as pd
import os
import networkx as nx
import scipy
import logging
import numpy as np
import argparse
import pickle
import itertools

from typing import Dict, Union, List, Tuple, Set
from statsmodels.stats.multitest import fdrcorrection

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

# digest Apocrita job name to give a log file
parser = argparse.ArgumentParser()
parser.add_argument('--gwas_disease', type=str)
parser.add_argument('--proteomics_disease', type=str)
parser.add_argument('--proteomics_algo', type=str, choices=['firth', 'cox'])
parser.add_argument('--tissue_index', type=int, default=1)  # SLURM_TASK_ID

args = parser.parse_args()  # parse arguments
tissue_index = args.tissue_index - 1
GWAS_DISEASE = args.gwas_disease
PROTEOMICS_DISEASE = args.proteomics_disease
ALGO = args.proteomics_algo

# select tissue
tissues = sorted([x.name for x in os.scandir(f"~/data/internal/knowledge_graph/tissue_graphs") if x.is_dir()])
TISSUE = tissues[tissue_index]

JOB_NAME = f"gwas{GWAS_DISEASE}_{ALGO}{PROTEOMICS_DISEASE}"
LOG_PATH = ''
LOGGING_PATH = f"{LOG_PATH}/gene_protein_links_{JOB_NAME}_{TISSUE}.log"

# In a module, __name__ = module's name in the Python package namespace
# __name__ is in the 'name' argument. logger with name foo is the parent of a logger with name foo.bar etc.
logger = logging.getLogger(f"GWAS-proteomics links for {TISSUE}: {JOB_NAME}")
logger.setLevel(logging.INFO)
logger.propagate = False

logger.handlers.clear()  # clear existing handlers, including streaming to console

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


def get_proteomics_genes(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> Set[str]:

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

    return genes


def get_gwas_magma_genes(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> Set[str]:

    logger.info(f"Pulling genes associated to {disease_code} through GWAS and MAGMA.")

    genes = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # iterate through the edges (there may be multiple edges/edge types)
            for edge_data in edges_data.values():
                # select edges which represent a MAGMA link
                if edge_data.get('edge_type') == 'magma_associated_with':
                    genes.add(neighbour)

    logger.info(f"Pulling GWAS MAGMA genes complete. Pulled {len(genes)} genes.")

    return genes


def get_gwas_eqtl_genes(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> Set[str]:
    logger.info(f"Pulling genes associated to {disease_code} through GWAS and eQTLs.")

    genes = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # iterate through the edges (there may be multiple edges/edge types)
            for edge_data in edges_data.values():
                # select edges which represent a GWAS link
                if edge_data.get('edge_type') == 'gwas_associated_with':
                    variant = neighbour
                    # now check genes linked to this variant - there can be many genes linked to the variant
                    for gene in full_graph.neighbors(variant):
                        variant_gene_edges = full_graph.get_edge_data(variant, gene)
                        if variant_gene_edges:
                            for variant_gene_edge in variant_gene_edges.values():
                                if variant_gene_edge.get('edge_type') == 'eqtl_for':
                                    genes.add(gene)

    logger.info(f"Pulling GWAS eQTL genes complete. Pulled {len(genes)} genes.")

    return genes


def get_gwas_pqtl_genes(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> Set[str]:
    logger.info(f"Pulling genes associated to {disease_code} through GWAS and pQTLs.")

    genes = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # iterate through the edges (there may be multiple edges/edge types)
            for edge_data in edges_data.values():
                # select edges which represent a GWAS link
                if edge_data.get('edge_type') == 'gwas_associated_with':
                    variant = neighbour
                    # now check genes linked to this variant - there can be many genes linked to the variant
                    for gene in full_graph.neighbors(variant):
                        variant_gene_edges = full_graph.get_edge_data(variant, gene)
                        if variant_gene_edges:
                            for variant_gene_edge in variant_gene_edges.values():
                                if variant_gene_edge.get('edge_type') == 'pqtl_for':
                                    genes.add(gene)

    logger.info(f"Pulling GWAS pQTL genes complete. Pulled {len(genes)} genes.")

    return genes


def get_gwas_variants(full_graph: nx.MultiDiGraph, disease_code: str = 'E11') -> Set[str]:
    logger.info(f"Pulling genes associated to {disease_code} through GWAS and pQTLs.")

    variants = set()
    # for each node neighbouring the given disease:
    for neighbour in full_graph.predecessors(disease_code):
        # get the edges connecting the neighbouring node and the disease
        edges_data = full_graph.get_edge_data(neighbour, disease_code)
        # if there are edges:
        if edges_data:
            # iterate through the edges (there may be multiple edges/edge types)
            for edge_data in edges_data.values():
                # select edges which represent a GWAS link
                if edge_data.get('edge_type') == 'gwas_associated_with':
                    variants.add(neighbour)

    logger.info(f"Pulling GWAS variants complete. Pulled {len(variants)} variants.")

    return variants


########################################################################################################################
#
# PULL NODE ATTRIBUTES
#
########################################################################################################################


def pull_secretome_information(full_graph: nx.MultiDiGraph, proteomics_genes: Set[str]) -> Dict[str, Set[str]]:
    # secreted at all, secreted to blood, not secreted
    secreted = {n for n, d in full_graph.nodes(data='feature_secretome_location')
                if isinstance(d, str)} & proteomics_genes
    secreted_to_blood = {n for n, d in full_graph.nodes(data='feature_secretome_location')
                         if 'blood' in str(d)} & proteomics_genes
    non_secreted = {n for n, d in full_graph.nodes(data='feature_secretome_location')
                    if d is None} & proteomics_genes

    return {'secreted': secreted, 'secreted_to_blood': secreted_to_blood, 'not_secreted': non_secreted}


########################################################################################################################
#
# ENRICHMENT FUNCTIONS
#
########################################################################################################################


def overrepresentation(test_set, background_annotations, entity_col_id, annotation_col_id):
    """
    Fisher's exact test for over-representation enrichment analysis.
    :param test_set: list of entity IDs to test for annotation enrichment
    :param background_annotations: DataFrame linking annotations to all entities in the background set. Must contain <entity_col_id> and <annotation_col_id>
    :param entity_col_id: column name in background_annotations for entity IDs, e.g. UniProt IDs for proteins
    :param annotation_col_id: column name in background_annotations for annotation IDs, e.g. R-HSA-xx for Reactome pathways
    :return: DataFrame listing enrichment p-values for each annotation present in the test_set
    """

    background_annotations = background_annotations.drop_duplicates([entity_col_id, annotation_col_id],
                                                                    ignore_index=True)

    n_total = len(background_annotations[entity_col_id].unique())

    # ensure test set is a strict subset of annotated background entities
    # otherwise entities in the test set could be not in the background
    test_annotations = background_annotations[background_annotations[entity_col_id].isin(test_set)].copy()
    # number of test entities that have at least one annotation
    n_test = len(test_annotations[entity_col_id].unique())
    # number of entities from background per annotation
    n_total_annotated_series = background_annotations[annotation_col_id].value_counts()
    # number of test entities per annotation
    n_test_annotated_series = test_annotations[annotation_col_id].value_counts()

    results = pd.DataFrame()
    for annotation in test_annotations[annotation_col_id].unique():
        n_total_annotated = n_total_annotated_series[annotation]
        n_test_annotated = n_test_annotated_series[annotation]

        # set up contingency table [[a, b], [c, d]] where columns are in-test-set +/-, rows are annotation +/-
        #                    In Test Set
        #                -------------------
        #                |   +   |   -     |
        #  --------------|-------|---------|
        #  Annotation +  |   a   |   b     |  n_total_annotated
        #  Annotation -  |   c   |   d     |
        #  --------------|-------|---------|
        #                 n_test              n_total

        a = n_test_annotated
        b = n_total_annotated - n_test_annotated
        c = n_test - n_test_annotated
        d = n_total - n_test - b

        fisher_test = scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        results = pd.concat([results, pd.DataFrame([{'annotation': annotation,
                                                     'odds_ratio': fisher_test.statistic,
                                                     'pval': fisher_test.pvalue,
                                                     'n_total_annotated': n_total_annotated,
                                                     'n_test_annotated': n_test_annotated}])])

    return results


def slim_annotations(annotations, annotation_col_id, entity_col_id,
                     lower_info_content_bound=None, upper_info_content_bound=None,
                     aspect_col_id=None):
    # set info_content bounds if not supplied
    # absolute thresholds taken from Donertas 2021 and converted to natural units of information by comparing against
    #  the number of genes annotated in the GO after filtering evidence codes and removing NOT relations
    if lower_info_content_bound is None:
        lower_info_content_bound = -np.log(500 / 19026)
    if upper_info_content_bound is None:
        upper_info_content_bound = -np.log(10 / 19026)

    # count annotation frequency, within aspect if specified#
    if aspect_col_id is not None:
        annotations_to_return = pd.DataFrame()
        for aspect in annotations[aspect_col_id].unique():
            aspect_df = annotations[annotations[aspect_col_id] == aspect].copy()
            # calcuate information content of all annotations
            annotation_info_content = -np.log(aspect_df[annotation_col_id].value_counts() /
                                              len(aspect_df[entity_col_id].unique()))
            chosen_annotations = (lower_info_content_bound < annotation_info_content) & (
                    annotation_info_content < upper_info_content_bound)
            aspect_annotations_to_return = aspect_df[
                aspect_df[annotation_col_id].isin(chosen_annotations[chosen_annotations].index)]
            annotations_to_return = pd.concat([annotations_to_return, aspect_annotations_to_return])
    else:
        # calcuate information content of all annotations
        annotation_info_content = -np.log(annotations[annotation_col_id].value_counts() /
                                          len(annotations[entity_col_id].unique()))

        chosen_annotations = (lower_info_content_bound < annotation_info_content) & (
                annotation_info_content < upper_info_content_bound)

        annotations_to_return = annotations[
            annotations[annotation_col_id].isin(chosen_annotations[chosen_annotations].index)]

    return annotations_to_return


def prep_annotation_files(raw_annotation_file):
    prepped_annotation_file = raw_annotation_file[
        (raw_annotation_file[6] == 'TAS') & (raw_annotation_file[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                             axis=1).set_axis(
        ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()

    return prepped_annotation_file


########################################################################################################################
#
# SIGNIFICANCE LABELLING
#
########################################################################################################################


def significance_labelling(input_data, pvalue_col='pval', alpha=0.05, fdr_permissive_group_by=None,
                           pair_mirroring=False,
                           pair_columns=None):
    input_data = input_data.copy()

    # if the input_data contains both permutations of a pair whose data are the same, i.e. not independent tests,
    #  then remove the duplication to accurately capture the number of independent tests
    if not pair_mirroring:
        processing_data = input_data.copy()
    else:
        input_data['pair_sorted'] = input_data.apply(
            lambda row: '_'.join(sorted([row[pair_columns[0]], row[pair_columns[1]]])), axis=1)
        processing_data = input_data.drop_duplicates(['pair_sorted', pvalue_col], ignore_index=True).copy()

    # add significance labels
    processing_data.loc[processing_data[pvalue_col] < alpha, 'nom_sig'] = 1
    processing_data['fdr_sig'], processing_data['pval_fdr_corrected'] = fdrcorrection(processing_data[pvalue_col],
                                                                                      alpha=alpha)
    processing_data.loc[processing_data[pvalue_col] < alpha / len(processing_data), 'bonf_sig'] = 1

    # the more generous FDR verion (correcting within groups only, rather than across the whole dataset)
    if fdr_permissive_group_by is not None:
        processing_data_processed = pd.DataFrame()
        for group in processing_data[fdr_permissive_group_by].unique():
            cluster_results = processing_data[processing_data[fdr_permissive_group_by] == group].copy()
            cluster_results['fdr_sig_permissive'], cluster_results['pval_fdr_corrected_permissive'] = fdrcorrection(
                cluster_results[pvalue_col], alpha=alpha)
            processing_data_processed = pd.concat([processing_data_processed, cluster_results])
    else:
        processing_data_processed = processing_data.copy()
    processing_data_processed = processing_data_processed.replace(True, 1).replace(False, 0).fillna(0)

    processing_data_processed['nom_sig'] = processing_data_processed['nom_sig'].astype(int)
    processing_data_processed['bonf_sig'] = processing_data_processed['bonf_sig'].astype(int)

    # if the input_data contains both permutations of a pair whose data are the same, i.e. not independent tests,
    #  restore both permutations of the pair to the output_data
    if not pair_mirroring:
        output_data = processing_data_processed.copy()
    else:
        output_data_med = pd.concat([processing_data_processed,
                                     processing_data_processed.assign(disease1=processing_data_processed['disease2'],
                                                                      disease2=processing_data_processed['disease1'])])
        output_data = pd.DataFrame()
        for pair_sorted in output_data_med['pair_sorted'].unique():
            output_data = pd.concat([output_data, output_data_med[output_data_med['pair_sorted'] == pair_sorted]])
        output_data.drop(columns='pair_sorted', inplace=True)
        # output_data = processing_data_processed.drop(columns=pair_columns).merge(input_data, on=[x for x in input_data.columns if x not in pair_columns], how='right').drop(columns='pair_sorted')

    return output_data


########################################################################################################################
#
# FUNCTIONS FOR SIMILARITY MEASURES
#
########################################################################################################################


# direct overlap - do the GWAS and proteomics sets share genes, and to significance?
def direct_overlap(
        full_graph: nx.MultiDiGraph, gwas_genes: Set[str], proteomics_genes: Set[str]) \
        -> Dict[str, Union[int, float, Set[str]]]:
    logger.info(f"Evaluating direct overlap between GWAS and proteomic gene sets.")

    all_genes = {n for n, d in full_graph.nodes(data='type') if d == 'gene'}
    both = gwas_genes & proteomics_genes
    only_gwas = gwas_genes - proteomics_genes
    only_proteomics = proteomics_genes - gwas_genes

    # set up contingency table [[a, b], [c, d]]
    #                    Proteomics
    #                -------------------
    #                |   +   |   -     |
    #  --------------|-------|---------|
    #        GWAS +  |   a   |   b     |  len(gwas_genes)
    #        GWAS -  |   c   |   d     |
    #  --------------|-------|---------|
    #           len(proteomics_genes)     n_total

    a = len(both)
    b = len(only_gwas)
    c = len(only_proteomics)
    d = len(all_genes) - len(proteomics_genes) - b

    fisher_test = scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')

    return_dict = {'fisher_greater_pval': fisher_test.pvalue,
                   'shared_genes': both,
                   'n_shared_genes': a}

    logger.info(f"Evaluation of direct overlap complete: {a} genes shared with greater_pval {fisher_test.pvalue}.")

    return return_dict


def get_shortest_path(full_graph: nx.MultiDiGraph, source: str, target: str) -> List[str]:
    # get the shortest path between two nodes, represented as [source, edge, node1, edge, ..., target]
    # nx.shortest_path returns [source, node1, node2, ..., target]
    if not nx.has_path(full_graph, source=source, target=target):
        return []
    shortest_path_nodes = nx.shortest_path(full_graph, source=source, target=target)
    shortest_path_nodes_edges = [f"node:{shortest_path_nodes[0]}"]
    for i in range(len(shortest_path_nodes) - 1):
        single_hop_edges = full_graph.get_edge_data(shortest_path_nodes[i], shortest_path_nodes[i + 1])
        for edge_data in single_hop_edges.values():
            shortest_path_nodes_edges.append(f"edge:{edge_data.get('edge_type')}")
        shortest_path_nodes_edges.append(f"node:{shortest_path_nodes[i + 1]}")

    return shortest_path_nodes_edges


# shortest links that exist between the two sets (exclude genes in both sets)
def shortest_paths_between_gwas_proteomics_sets(
        full_graph: nx.MultiDiGraph, gwas_nodes: Set[str], proteomics_genes: Set[str]) -> Dict[str, List[str]]:
    logger.info(f"Finding shortest paths between each GWAS gene and each proteomic gene.")
    only_gwas = gwas_nodes - proteomics_genes
    only_proteomics = proteomics_genes - gwas_nodes
    gene_pairs = list(itertools.product(list(only_gwas), list(only_proteomics)))
    gene_pairs_shortest_paths = {
        gene_pair: get_shortest_path(full_graph=full_graph, source=gene_pair[0], target=gene_pair[1])
        for gene_pair in gene_pairs}

    # remove empty paths, i.e. no path exists between source and target
    gene_pairs_shortest_paths = {k: v for k, v in gene_pairs_shortest_paths.items() if len(v) > 0}

    # order by len(value lists) - i.e. shortest paths first, in ascending order
    gene_pairs_shortest_paths = dict(sorted(gene_pairs_shortest_paths.items(), key=lambda item: len(item[1])))

    if len(gene_pairs_shortest_paths) == 0:
        logger.info(f"No shortest paths found.")
        return {}

    logger.info(f"Finding {len(gene_pairs_shortest_paths)} shortest paths complete. "
                f"Data returned from shortest path "
                f"({len([x for x in gene_pairs_shortest_paths[list(gene_pairs_shortest_paths)[0]] if 'node' in x]) - 1} "
                f"hops) ascending.")

    return gene_pairs_shortest_paths


def execution():
    # TISSUE = 'Adipose_Subcutaneous'
    # GWAS_DISEASE = 'G30'
    # PROTEOMICS_DISEASE = 'G30'
    # disease = 'G20'

    tissue_graph = load_graph(tissue=TISSUE)
    # tissue_graph = nx.read_graphml('ukbiobank/cellular_network/tissue_graphs/Adipose_Subcutaneous/graph.graphml')
    gene_nodes = sorted([n for n, d in list(tissue_graph.nodes(data='type')) if d == 'gene'])
    gene_only_graph = tissue_graph.subgraph(gene_nodes).copy()

    # load proteomics set - can use multiple proteomics diseases who then have their assoc proteins combined
    proteomic_genes = set()
    for proteomics_disease in PROTEOMICS_DISEASE.split('-'):
        proteomic_genes |= get_proteomics_genes(full_graph=tissue_graph, disease_code=proteomics_disease)
    proteomic_secretome_info = pull_secretome_information(full_graph=tissue_graph, proteomics_genes=proteomic_genes)

    # load GWAS gene sets
    gwas_eqtl_genes = set()
    gwas_pqtl_genes = set()
    gwas_magma_genes = set()
    gwas_variants = set()
    for gwas_disease in GWAS_DISEASE.split('-'):
        gwas_eqtl_genes |= get_gwas_eqtl_genes(full_graph=tissue_graph, disease_code=gwas_disease)
        gwas_pqtl_genes |= get_gwas_pqtl_genes(full_graph=tissue_graph, disease_code=gwas_disease)
        gwas_magma_genes |= get_gwas_magma_genes(full_graph=tissue_graph, disease_code=gwas_disease)
        gwas_variants |= get_gwas_variants(full_graph=tissue_graph, disease_code=gwas_disease)
    gwas_gene_sets = {'eqtl': gwas_eqtl_genes, 'pqtl': gwas_pqtl_genes, 'magma': gwas_magma_genes,
                      'all_gwas': gwas_eqtl_genes | gwas_pqtl_genes | gwas_magma_genes}

    # check gene sets have stuff
    if len(proteomic_genes) == 0:
        logger.info(f"There are no proteomics ({ALGO}) genes for {PROTEOMICS_DISEASE}.")
        raise ValueError(f"Proteomics genes are required to run this pipeline.")
    if len(gwas_eqtl_genes | gwas_pqtl_genes | gwas_magma_genes) == 0:
        logger.info(f"There are no GWAS genes for {GWAS_DISEASE}.")
        raise ValueError(f"GWAS genes are required to run this pipeline.")

    # # import gene embeddings
    # logger.info(f"Importing gene embeddings for {TISSUE}.")
    # with open(f"{SCRATCH_APOCRITA_DIR}/tissue_graphs/{TISSUE}/n2v_data_interacts-with_coexpresses-with.pkl",
    #           'rb') as pickle_file:
    #     embeddings = pickle.load(pickle_file)
    # logger.info(f"Importing embeddings complete.")

    export_dict = {}
    for gwas_gene_set_name, gwas_gene_set in gwas_gene_sets.items():

        logger.info(f"------Executing analyses for {gwas_gene_set_name} GWAS gene set.------")

        if len(gwas_gene_set) == 0:
            logger.info(f"There are no GWAS ({gwas_gene_set_name}) for {GWAS_DISEASE}. Skipping to next GWAS gene set.")
            continue

        # direct overlap
        proteomic_gwas_direct_overlap = direct_overlap(full_graph=tissue_graph,
                                                       gwas_genes=gwas_gene_set,
                                                       proteomics_genes=proteomic_genes)

        # shortest paths between sets
        shortest_paths_gwas_genes_proteomic_genes = shortest_paths_between_gwas_proteomics_sets(
            full_graph=gene_only_graph,
            gwas_nodes=gwas_gene_set,
            proteomics_genes=proteomic_genes)

        # make results dict
        export_dict[gwas_gene_set_name] = {
            'direct_overlap': proteomic_gwas_direct_overlap,
            'shortest_paths_between_sets': shortest_paths_gwas_genes_proteomic_genes,
            'gwas_set': gwas_gene_set,
            'proteomic_set': proteomic_genes
        }

    # export all results - pickle
    logger.info(f"Exporting results to '~/data/internal/knowledge_graph/gene_protein_graph_analysis/{JOB_NAME}/graph_analysis_pickles/{TISSUE}.pkl'.")
    os.makedirs(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{JOB_NAME}/graph_analysis_pickles", exist_ok=True)
    with open(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{JOB_NAME}/graph_analysis_pickles/{TISSUE}.pkl", 'wb') as pickle_file:
        pickle.dump(export_dict, pickle_file)

    logger.info(f"Export and execution complete.")


if __name__ == '__main__':
    execution()
