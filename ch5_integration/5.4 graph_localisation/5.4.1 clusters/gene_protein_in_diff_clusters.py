"""
Detect clusters in the graph, and test whether GWAS genes and proteomics genes are enriched in each of the clusters.
Then, test whether there is more co-occurrence of GWAS and proteomic enrichment in the same clusters than expected,
 for each tissue-disease combination.
Overall, tests the hypothesis that GWAS gene and proteomics genes are separated in the cellular network.
"""
import pandas as pd
import os
import networkx as nx
import scipy
import logging
import numpy as np
import statsmodels.stats.multitest
import argparse

from typing import Dict, List, Tuple, Set

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

# digest Apocrita job name to give a log file
parser = argparse.ArgumentParser()
# parser.add_argument('--jobname', type=str, required=True)
parser.add_argument('--gene_edges_to_use', type=str, choices=['interacts_with', 'coexpresses_with', 'both'],
                    default='both')
parser.add_argument('--proteomics_algo_to_use', type=str, choices=['firth', 'cox'], default='firth')
args = parser.parse_args()  # parse arguments
# job_name = args.jobname
gene_edges_to_use = args.gene_edges_to_use
ALGO = args.proteomics_algo_to_use

if gene_edges_to_use == 'both':
    gene_edges_used = ['interacts_with', 'coexpresses_with']
else:
    gene_edges_used = [gene_edges_to_use]
job_name = '_'.join([x.replace('_', '-') for x in gene_edges_used])
# add which proteomics links are being used: firth or cox
job_name = f"{job_name}_{ALGO}"
LOG_PATH = ''
LOGGING_PATH = f"{LOG_PATH}/gene_protein_in_diff_clusters_{job_name}.log"

# In a module, __name__ = module's name in the Python package namespace
# __name__ is in the 'name' argument. logger with name foo is the parent of a logger with name foo.bar etc.
# logger = logging.getLogger('GWAS vs proteomic genes in the cellular network')
logger = logging.getLogger(job_name)
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

    graph = nx.read_graphml(f"~/data/internal//tissue_graphs/{tissue}/graph.graphml")

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
                    f"There are {len(edge_data)} edges between node {neighbour} and disease {disease_code}.")
            # iterate through the edges
            for edge_data in edges_data.values():
                # select edges which represent a proteomics Firth/Cox link
                if edge_data.get('edge_type') == f"proteomics_{ALGO}_associated_with":
                    genes.add(neighbour)

    logger.info(f"Pulling proteomics genes complete.")

    return list(genes)


def get_gwas_genes(full_graph: nx.MultiDiGraph, gene_gwas_relationship: List[str], disease_code: str = 'E11') -> List[str]:

    logger.info(f"Pulling genes associated to {disease_code} through GWAS (including MAGMA).")

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

    logger.info('Pulling GWAS genes complete.')

    return list(genes)


########################################################################################################################
#
# COMMUNITY DETECTION: detection, enrichment
#
########################################################################################################################


def create_homo_gene_subgraph(full_graph: nx.MultiDiGraph, inter_gene_edges: List[str]) -> nx.Graph:
    logger.info(f"Creating a homogeneous gene subgraph with edge(s): {' and '.join(inter_gene_edges)}.")

    # converting to sorted list to ensure ordering of the nodes in the subgraph is consistent between runs
    genes = sorted(list({n for n, d in full_graph.nodes(data=True) if d.get("type") == 'gene'}))
    selected_edges = sorted([(u, v) for u, v, k, d in full_graph.edges(keys=True, data=True)
                             if u in genes and v in genes and d.get("edge_type") in inter_gene_edges])
    subgraph = nx.Graph()
    subgraph.add_nodes_from(genes)
    subgraph.add_edges_from(selected_edges)

    logger.info('Homogeneous gene subgraph created. Converting node labels to integer IDs.')

    # ordering argument forces a consistent mapping between the original node_ids and the integer IDs
    subgraph_int_ids = nx.convert_node_labels_to_integers(subgraph, label_attribute='original_id', ordering='sorted')

    logger.info('Conversion complete.')

    return subgraph_int_ids


def louvain_clustering(full_graph, inter_gene_edges: List[str]) -> List[set]:
    logger.info('Louvain clustering of the homogeneous gene subgraph.')

    inter_gene_edge_check = set(inter_gene_edges) - {'interacts_with', 'coexpresses_with'}
    if len(inter_gene_edge_check) > 0:
        logger.error(f"Invalid gene-gene edges: {inter_gene_edge_check}. "
                     f"Must be 'interacts_with' and/or 'coexpresses_with'.")
        raise ValueError(f"Invalid gene-gene edges: {inter_gene_edge_check}. "
                         f"Must be 'interacts_with' and/or 'coexpresses_with'.")

    louvain_subgraph = create_homo_gene_subgraph(full_graph=full_graph, inter_gene_edges=inter_gene_edges)
    louvain_clusters = nx.community.louvain_communities(louvain_subgraph, weight=None, seed=42)

    # mapping from integer IDs to original IDs
    int_to_og = nx.get_node_attributes(louvain_subgraph, 'original_id')

    # map clusters back to original IDs
    louvain_clusters_orig = [set(int_to_og[i] for i in cluster) for cluster in louvain_clusters]

    logger.info('Clustering complete.')

    return louvain_clusters_orig


def fisher_enrichment(
        axis1_gene_set: Set[str], axis2_gene_set: Set[str], background_gene_set: Set[str]) -> Tuple[float, float]:
    # set up contingency table [[a, b], [c, d]] where columns are in-test-set +/-, rows are annotation +/-
    #                  GWAS/proteomics
    #                -------------------
    #                |   +   |   -     |
    #  --------------|-------|---------|
    #     Cluster +  |   a   |   b     |
    #     Cluster -  |   c   |   d     |
    #  --------------|-------|---------|
    #

    # axis1 = GWAS/proteomics
    # axis2 = cluster
    # background = all genes

    a = len(axis1_gene_set & axis2_gene_set)
    b = len(axis2_gene_set - axis1_gene_set)
    c = len(axis1_gene_set - axis2_gene_set)
    d = len(background_gene_set) - a - b - c

    fisher_test = scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')

    # if a row or column has only zeros, the table becomes degenerate and the odds ratio is set to np.nan
    # pval may still be informative
    if (a + b) * (c + d) * (a + c) * (b + d) == 0:
        return np.nan, fisher_test.pvalue

    return fisher_test.statistic, fisher_test.pvalue


def cluster_enrichment(clusters: List[Set[str]], gwas_genes: List[str], proteomics_genes: List[str]) -> pd.DataFrame:
    gwas_genes_set = set(gwas_genes)
    proteomics_genes_set = set(proteomics_genes)

    # get all genes (assumed that all genes have ended up assigned to a cluster)
    all_genes = set.union(*clusters)

    logger.info('Testing clusters for enrichment in genes associated to disease by GWAS and proteomics.')

    cluster_results = pd.DataFrame()
    for cluster_count, cluster in enumerate(clusters):

        # if there are fewer than 3 nodes in a cluster, that doesn't really count as a cluster
        #  and should not be tested for enrichment
        if len(cluster) < 3:
            continue

        gwas_or, gwas_pval = fisher_enrichment(axis1_gene_set=gwas_genes_set, axis2_gene_set=cluster,
                                               background_gene_set=all_genes)
        proteomics_or, proteomics_pval = fisher_enrichment(axis1_gene_set=proteomics_genes_set, axis2_gene_set=cluster,
                                                           background_gene_set=all_genes)

        cluster_results = pd.concat([cluster_results,
                                     pd.DataFrame({'cluster': [cluster_count] * 2,
                                                   'n_cluster_genes': len(cluster),
                                                   'or': [gwas_or, proteomics_or],
                                                   'pval': [gwas_pval, proteomics_pval],
                                                   'gene_set': ['gwas', 'proteomics']})])

    # remove non-tests
    cluster_results = cluster_results.dropna(subset='pval', axis=0, ignore_index=True)

    logger.info('Enrichment testing complete. Significance labelling.')

    # label with nominal significance
    cluster_results['nom_sig'] = (cluster_results['pval'] < 0.05).astype(int)

    # FDR correct
    cluster_results['fdr_sig'], cluster_results['pval_fdr_corrected'] = statsmodels.stats.multitest.fdrcorrection(
        cluster_results['pval'], alpha=0.05)
    cluster_results['fdr_sig'] = cluster_results['fdr_sig'].astype(int)

    # pivot wider to be: cluster, and gwas/proteomic odds ratio, pval, fdr_sig, pval_fdr_corrected
    cluster_results_wide = (cluster_results[cluster_results['gene_set'] == 'gwas']
                            .rename(columns={'or': 'gwas_or', 'pval': 'gwas_pval', 'nom_sig': 'gwas_nom_sig',
                                             'fdr_sig': 'gwas_fdr_sig',
                                             'pval_fdr_corrected': 'gwas_pval_fdr_corrected'})
                            .drop(columns='gene_set')
                            .merge(cluster_results[cluster_results['gene_set'] == 'proteomics']
                                   .rename(columns={'or': 'proteomics_or', 'pval': 'proteomics_pval',
                                                    'nom_sig': 'proteomics_nom_sig', 'fdr_sig': 'proteomics_fdr_sig',
                                                    'pval_fdr_corrected': 'proteomics_pval_fdr_corrected'})
                                   .drop(columns='gene_set'), on=['cluster', 'n_cluster_genes'], how='outer'))

    return cluster_results_wide


def statistically_evaluate_clusters(label_counts_dict: Dict[str, pd.Series]) -> pd.DataFrame:
    summary_df = pd.DataFrame()
    for sig_level, label_counts in label_counts_dict.items():
        # see if there is depletion for gwas and proteomics proteins existing in the same clusters
        # if so, this evidences the hypothesis of gwas and proteomics separation in the graph
        contingency_table = [[label_counts.get('both', 0), label_counts.get('gwas', 0)],
                             [label_counts.get('proteomics', 0), label_counts.get('neither', 0)]]
        fisher_less_pval = scipy.stats.fisher_exact(contingency_table, alternative='less').pvalue
        fisher_greater_pval = scipy.stats.fisher_exact(contingency_table, alternative='greater').pvalue
        summary_df = pd.concat([
            summary_df,
            pd.DataFrame([{'sig_level': sig_level, 'fisher_less_pval': fisher_less_pval,
                           'fisher_greater_pval': fisher_greater_pval} |
                          {k + '_count': v for k, v in label_counts.to_dict().items()}])])

    return summary_df  # sig_level, fisher_greater_pval, fisher_less_pval, x_count (x = gwas, proteomics, both, neither)


def evaluate_clusters(cluster_results_wide: pd.DataFrame):
    logger.info('Evaluating whether any clusters share GWAS and proteomics genes.')

    # label rows with gene/protein set specificity
    cluster_results_wide['fdr_specificity_label'] = 'neither'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_fdr_sig'] > cluster_results_wide[
            'proteomics_fdr_sig'], 'fdr_specificity_label'] = 'gwas'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_fdr_sig'] < cluster_results_wide[
            'proteomics_fdr_sig'], 'fdr_specificity_label'] = 'proteomics'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_fdr_sig'] * cluster_results_wide[
            'proteomics_fdr_sig'] == 1, 'fdr_specificity_label'] = 'both'

    label_counts_fdr = cluster_results_wide['fdr_specificity_label'].value_counts()

    logger.info(f"Out of {len(cluster_results_wide)} clusters: \n"
                f"- {label_counts_fdr.get('gwas', 0)} solely GWAS-enriched (FDR < 0.05) clusters \n"
                f"- {label_counts_fdr.get('proteomics', 0)} solely proteomics-enriched (FDR < 0.05) clusters \n"
                f"- {label_counts_fdr.get('both', 0)} GWAS- and proteomics-enriched (FDR < 0.05) clusters.")

    # label rows with gene/protein set specificity
    cluster_results_wide['specificity_label'] = 'neither'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_nom_sig'] > cluster_results_wide['proteomics_nom_sig'], 'specificity_label'] = 'gwas'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_nom_sig'] < cluster_results_wide[
            'proteomics_nom_sig'], 'specificity_label'] = 'proteomics'
    cluster_results_wide.loc[
        cluster_results_wide['gwas_nom_sig'] * cluster_results_wide[
            'proteomics_nom_sig'] == 1, 'specificity_label'] = 'both'

    label_counts = cluster_results_wide['specificity_label'].value_counts()

    logger.info(f"Out of {len(cluster_results_wide)} clusters: \n"
                f"- {label_counts.get('gwas', 0)} solely GWAS-enriched (p < 0.05) clusters \n"
                f"- {label_counts.get('proteomics', 0)} solely proteomics-enriched (p < 0.05) clusters \n"
                f"- {label_counts.get('both', 0)} GWAS- and proteomics-enriched (p < 0.05) clusters.")

    logger.info('Testing whether there is statistically significant depletion/enrichment of GWAS- and '
                'proteomics-enriched clusters.')
    # summary_data: sig_level, fisher_greater_pval, fisher_less_pval, x_count (x = gwas, proteomics, both, neither)
    summary_data = statistically_evaluate_clusters({'fdr': label_counts_fdr, 'nominal': label_counts})

    logger.info('Evaluation complete.')

    return cluster_results_wide, summary_data


def fdr_correct_collated_summaries(collated_summaries: pd.DataFrame) -> pd.DataFrame:

    def fdr_sub_func(grouped_df):

        grouped_df_working = grouped_df.copy()
        grouped_df_working['fisher_less_fdr_sig'], _ = statsmodels.stats.multitest.fdrcorrection(
            grouped_df_working['fisher_less_pval'], alpha=0.05)
        grouped_df_working['fisher_greater_fdr_sig'], _ = statsmodels.stats.multitest.fdrcorrection(
            grouped_df_working['fisher_greater_pval'], alpha=0.05)

        grouped_df_working['fisher_less_fdr_sig'] = grouped_df_working['fisher_less_fdr_sig'].astype(int)
        grouped_df_working['fisher_greater_fdr_sig'] = grouped_df_working['fisher_greater_fdr_sig'].astype(int)

        return grouped_df_working

    fdr_corrected = collated_summaries.groupby('sig_level').apply(fdr_sub_func).reset_index(drop=True)

    return fdr_corrected


########################################################################################################################
#
# ADMIN FUNCTIONS: prep clusters for export, label results for collation
#
########################################################################################################################


def prep_clusters_for_export(clusters: List[Set[str]]) -> pd.DataFrame:
    cluster_df = pd.DataFrame({
        'cluster': range(len(clusters)),
        'cluster_genes': ['|'.join(sorted(s)) for s in clusters]
    })

    return cluster_df


def label_results_for_collation(results_df: pd.DataFrame, cols_to_add: Dict[str, str]) -> pd.DataFrame:

    results_df_labelled = results_df.copy()
    for col_name, col_content in cols_to_add.items():
        results_df_labelled[col_name] = col_content

    return results_df_labelled


########################################################################################################################
#
# MAIN EXECUTIONS
#
########################################################################################################################


def execute_all_tissue_disease(edges_to_use: str) -> None:

    tissues = [x.name for x in os.scandir(f"~/data/internal/knowledge_graph/tissue_graphs") if x.is_dir()]
    disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    diseases = disease_info['icd10_three_letter'].to_list()

    if edges_to_use == 'both':
        edges_used = ['interacts_with', 'coexpresses_with']
    else:
        edges_used = [edges_to_use]

    gene_gwas_relationships_list = [['variant-linked'], ['magma-linked'], ['variant-linked', 'magma-linked']]

    full_enrichment_results = pd.DataFrame()
    full_enrichment_summaries = pd.DataFrame()
    for tissue in tissues:

        tissue_graph = load_graph(tissue=tissue)
        # subgraph + cluster
        clusters = louvain_clustering(full_graph=tissue_graph,
                                      inter_gene_edges=edges_used)
        cluster_df = prep_clusters_for_export(clusters=clusters)
        cluster_df.to_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/clusters_{job_name}.csv", index=False)

        for disease in diseases:

            for gene_gwas_relationship in gene_gwas_relationships_list:

                proteomics_genes = get_proteomics_genes(full_graph=tissue_graph, disease_code=disease)
                gwas_genes = get_gwas_genes(full_graph=tissue_graph, gene_gwas_relationship=gene_gwas_relationship,
                                            disease_code=disease)

                # enrichment
                cluster_enrichment_results = cluster_enrichment(clusters=clusters, gwas_genes=gwas_genes,
                                                                proteomics_genes=proteomics_genes)

                # evaluate enrichment
                enrichment_results_labelled, enrichment_results_summary = evaluate_clusters(cluster_enrichment_results)

                # label results for collation
                collation_label_dict = {'tissue': tissue, 'disease': disease,
                                        'gwas_gene_disease_relationship': '_'.join(gene_gwas_relationship)}
                enrichment_results_labelled = label_results_for_collation(enrichment_results_labelled,
                                                                          cols_to_add=collation_label_dict)
                enrichment_results_summary = label_results_for_collation(enrichment_results_summary,
                                                                         cols_to_add=collation_label_dict)

                full_enrichment_results = pd.concat([full_enrichment_results, enrichment_results_labelled])
                full_enrichment_summaries = pd.concat([full_enrichment_summaries, enrichment_results_summary])

            # remove for the full run
            # break
        # break

    # note that this is FDR corrected within each tissue-disease-geneGWAS combination
    full_enrichment_results.to_csv(
        f"~/data/internal/knowledge_graph/gene_protein_distances/clustering/results_{job_name}.csv",
        index=False)

    # label collated summaries with significance and export
    # label with nominal significance
    full_enrichment_summaries['fisher_less_nom_sig'] = (
            full_enrichment_summaries['fisher_less_pval'] < 0.05).astype(int)
    full_enrichment_summaries['fisher_greater_nom_sig'] = (
            full_enrichment_summaries['fisher_greater_pval'] < 0.05).astype(int)
    # FDR correct across all tissue-disease combinations, but within each sig_label and each direction
    full_enrichment_summaries_fdr = fdr_correct_collated_summaries(full_enrichment_summaries)
    # fill NAs in x_count columns
    count_cols = [c for c in full_enrichment_summaries_fdr.columns if '_count' in c]
    full_enrichment_summaries_fdr[count_cols] = full_enrichment_summaries_fdr[count_cols].fillna(0)
    full_enrichment_summaries_fdr.to_csv(
        f"~/data/internal/knowledge_graph/gene_protein_distances/clustering/summaries_{job_name}.csv",
        index=False)


if __name__ == '__main__':
    execute_all_tissue_disease(edges_to_use=gene_edges_to_use)
