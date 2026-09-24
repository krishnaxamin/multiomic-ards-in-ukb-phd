""" Collate numbers of nodes and edges, including of each type, for each graph. """
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.use('TkAgg')

tissues = [x.name for x in os.scandir('~/data/internal/knowledge_graph/tissue_graphs') if x.is_dir()]

""" Read in stats """
tissue_stats_list = []
for tissue in tissues:
    tissue_stats_list.append(
        pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/graph_stats.csv").assign(tissue=tissue))
tissue_stats = pd.concat(tissue_stats_list)

total_nodes_edges = tissue_stats.reset_index(drop=True).pivot(index='tissue', columns='stat', values='value').fillna(
    0)[['n_nodes', 'n_edges']].astype(int).reset_index()

""" Read in number of each node and edge type """
type_counts_list = []
for tissue in tissues:
    nodes = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/nodes.csv")
    edges = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/edges.csv")

    node_type_counts = nodes.value_counts('type').reset_index().assign(tissue=tissue)
    edge_type_counts = edges.value_counts('edge_type').reset_index().set_axis(['type', 'count'], axis=1).assign(
        tissue=tissue)

    type_counts_list.append(node_type_counts)
    type_counts_list.append(edge_type_counts)

type_counts = pd.concat(type_counts_list)
type_counts_wide = type_counts.reset_index(drop=True).pivot(index='tissue', columns='type', values='count').fillna(
    0).astype(int).reset_index()
type_counts_wide = type_counts_wide.merge(total_nodes_edges)
type_counts_wide['tissue'] = [
    f"{x.split('_')[0]} ({x.split('_')[1]} {' '.join([xx.lower() for xx in x.split('_')[2:]])})".replace(' )', ')')
    if '_' in x else x
    for x in type_counts_wide['tissue'].to_list()]
type_counts_wide['tissue'] = type_counts_wide['tissue'].replace('Small (Intestine terminal ileum)',
                                                                'Small Intestine (Terminal ileum)')
type_counts_wide = type_counts_wide[['tissue', 'n_nodes', 'disease', 'gene', 'variant', 'n_edges'] +
                                    sorted([x for x in type_counts_wide.columns if '_' in x])]
type_counts_wide.to_csv('~/data/internal/knowledge_graph/n_node_edge_types_per_tissue.csv', index=False)
