import pandas as pd
import pickle
import os
import collections
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot
import numpy as np

from typing import Dict, Any, Tuple, List
from utils.significance_labelling import significance_labelling
from utils.enrichment_analyses import prep_annotation_files, slim_annotations, overrepresentation, \
    reactome_parentage_analysis

""" Ensembl gene set """


def generate_ensembl_gene_set(current_version: str) -> pd.DataFrame:
    """
    Pulls the current Ensembl gene set: every gene, its coordinates and its ENSG ID in Ensembl.

    :param current_version: Current ensembl version, e.g. 115
    :type current_version: str
    :returns: Current Ensembl gene set
    :rtype: pd.DataFrame
    """
    ensembl_gene_set = pd.read_csv(
        'https://ftp.ensembl.org/pub/release-' + current_version + '/gtf/homo_sapiens/Homo_sapiens.GRCh38.' + current_version + '.gtf.gz',
        sep='\t', skiprows=5, header=None)
    ensembl_gene_set_only_genes = ensembl_gene_set[ensembl_gene_set[2] == 'gene']
    ensembl_gene_set_only_genes.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                                           'attribute']

    # process attribute column to generate multiple columns for each attribute listed
    regex_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    attributes = ensembl_gene_set_only_genes['attribute'].map(
        lambda s: dict(regex_pattern.findall(str(s).strip().strip('"').replace('""', '"')))
    ).apply(pd.Series)
    ensembl_gene_set_only_genes = pd.concat([ensembl_gene_set_only_genes, attributes], axis=1)

    # export
    return ensembl_gene_set_only_genes.drop(columns=['score', 'frame', 'attribute'])


""" Reactome enrichment """


def reactome_ora(reactome_annots: pd.DataFrame, tissue: str, test_gene_set: set[str],
                 background: None | list[str] = None) -> pd.DataFrame:
    # define background (all genes in tissue-specific graph)
    background_genes = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/nodes.csv")['node_id'].to_list()
    if background:
        background_genes = list(set(background_genes) & set(background))
    reactome_background = reactome_annots[reactome_annots['id'].isin(background_genes)].copy()

    # slim background annotation file (note results for this may differ between UniProt/NCBI/Ensembl mappings)
    reactome_background_slimmed = slim_annotations(annotations=reactome_background, entity_col_id='id',
                                                   annotation_col_id='reactome_annotation_id')

    ora_result = overrepresentation(test_set=list(test_gene_set),
                                    background_annotations=reactome_background_slimmed,
                                    entity_col_id='id',
                                    annotation_col_id='reactome_annotation_id')

    if len(ora_result) == 0:
        return pd.DataFrame()

    ora_result = significance_labelling(ora_result)
    ora_result_sig = ora_result[ora_result['fdr_sig'] == 1].copy()

    return ora_result_sig


""" UpSet scaling """


def scaled_upset_plot(upsetplot_data):
    # given upsetdata (generated from from_contents or similar), generate parameters that scale the upsetplot
    #  proportionally to visually fit everything properly

    # upsetplot_data will have index of form (True, True, False, ...)
    #  this gets the length of that boolean vector = number of sets there are
    n_sets = len(upsetplot_data.index.names)
    n_intersections = len(set(upsetplot_data.index))

    # figure size
    fig_width = np.clip(n_intersections * 0.6, 8, 40)
    fig_height = np.clip(n_sets * 0.5, 5, 20)

    # text size: shrinks as sets increase, floored at 5pt
    label_size = np.clip(12 - (n_sets - 5) * 0.3, 5, 12)

    # dot and line size: shrink as the matrix gets denser
    dot_size = np.clip(200 - (n_sets - 5) * 6, 20, 200)
    linewidth = np.clip(1.2 - (n_sets - 5) * 0.03, 0.4, 1.2)

    fig = plt.figure(figsize=(fig_width, fig_height))
    # totals_plot_elements=0 removes the horizontal per-set bar charts
    gene_set_dict_working_upsetplot = upsetplot.UpSet(upsetplot_data, element_size=None,
                                                      totals_plot_elements=0, show_counts=True)
    gene_set_dict_working_upsetplot.style_subsets(min_degree=n_sets,
                                                  facecolor='#1f77b4')  # colour by 1st tab10 colour

    upsetplot_axes_dict = gene_set_dict_working_upsetplot.plot(fig=fig)

    ax_matrix = upsetplot_axes_dict['matrix']
    for collection in ax_matrix.collections:
        if isinstance(collection, mpl.collections.PathCollection):
            collection.set_sizes([dot_size])
        elif isinstance(collection, mpl.collections.LineCollection):
            collection.set_linewidths(linewidth)
    # for collection in ax_matrix.collections:
    #     collection.set_sizes([dot_size])      # scatter circles
    # for line in ax_matrix.lines:
    #     line.set_linewidth(linewidth)         # connecting lines

    ax_matrix.tick_params(axis='y', labelsize=label_size)
    upsetplot_axes_dict['intersections'].tick_params(labelsize=label_size * 0.9)
    if upsetplot_axes_dict.get('totals') is not None:
        upsetplot_axes_dict['totals'].tick_params(labelsize=label_size * 0.9)

    fig.tight_layout()
    return fig, upsetplot_axes_dict


""" Specificity-related functions """


def plot_gene_set_specificity(all_test_results: Dict[str, Any], analysis: str):
    # plot (1) number of genes in all genes specific in any tissues or not specific in any tissues
    #      (2) for each gene set, whether any tissue-specific genes are shared between tissues (upset)

    specific_dict = collections.defaultdict(set)
    nonspecific_dict = collections.defaultdict(set)
    tissue_specific_genes_per_geneset_dict = collections.defaultdict(
        dict)  # {'pqtl': {tissue: set()}, 'proteomics': {tissue: set()}, ...}
    for tissue, tissue_dict in all_test_results.items():
        # read in relevant graph nodes info
        tissue_nodes = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/nodes.csv")
        specific_genes = set(tissue_nodes[tissue_nodes['feature_specificity_category'].notna()]['node_id'])
        for gene_set_name, gene_set_dict in tissue_dict.items():
            # collate genes across tissues in different gene sets that are tissue-specific or not tissue-specific
            specific_dict[gene_set_name] |= gene_set_dict['gwas_set'] & specific_genes
            nonspecific_dict[gene_set_name] |= gene_set_dict['gwas_set'] - specific_genes
            specific_dict['proteomics'] |= gene_set_dict['proteomic_set'] & specific_genes
            nonspecific_dict['proteomics'] |= gene_set_dict['proteomic_set'] - specific_genes

            # collate genes specific to each tissue in each gene set
            tissue_specific_genes_per_geneset_dict[gene_set_name][tissue] = list(
                gene_set_dict['gwas_set'] & specific_genes)
            tissue_specific_genes_per_geneset_dict['proteomics'][tissue] = list(gene_set_dict[
                                                                                    'proteomic_set'] & specific_genes)
    # collate number of genes per gene set that are specific to any tissue vs not specific to any tissue
    specific_nonspecific_df = pd.concat(
        [pd.DataFrame([{'gene_set': k, 'gene_count': len(v)}
                       for k, v in specific_dict.items()]).assign(category='Specific'),
         pd.DataFrame([{'gene_set': k, 'gene_count': len(v)}
                       for k, v in nonspecific_dict.items()]).assign(category='Not specific')
         ])
    gene_set_order = ['eqtl', 'pqtl', 'magma', 'all_gwas', 'proteomics']
    specific_nonspecific_df['gene_set'] = pd.Categorical(specific_nonspecific_df['gene_set'],
                                                         [x for x in gene_set_order
                                                          if x in specific_nonspecific_df['gene_set'].to_list()])
    # plot number of genes per gene set that are specific to any tissue vs not specific to any tissue
    plt.figure(figsize=(8, 6))
    sns.barplot(data=specific_nonspecific_df, x='gene_set', y='gene_count', hue='category',
                hue_order=['Specific', 'Not specific'], palette=(plt.cm.coolwarm(1.0), plt.cm.coolwarm(0.0)))
    plt.xlabel('Gene set')
    plt.ylabel('Gene count')
    plt.gca().yaxis.grid(True, which='major', color='grey', alpha=0.6)
    plt.gca().set_axisbelow(True)
    plt.legend(title='Specificity status')
    plt.tight_layout()
    plt.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/genes_in_gene_sets_specificity_status.svg")
    plt.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/genes_in_gene_sets_specificity_status.png")
    plt.close()

    # plot the intersections between tissue-specific genes in the different gene sets,
    #  i.e. intersection between liver- and lung-specific gene sets
    for gene_set, gene_set_dict in tissue_specific_genes_per_geneset_dict.items():

        # don't do the Upset if there's only one tissue with tissue-specific genes
        gene_set_dict_working = {k: v for k, v in gene_set_dict.items() if len(v) > 0}
        if len(gene_set_dict_working) == 1:
            single_tissue = list(gene_set_dict_working)[0]
            print(
                f"'{gene_set}' has {len(gene_set_dict_working[single_tissue])} genes specific to {single_tissue} only.")
            continue
        elif len(gene_set_dict_working) == 0:
            print(
                f"'{gene_set}' has no genes specific to any tissues.")
            continue
        gene_set_dict_working_upset = upsetplot.from_contents(gene_set_dict_working)

        print(f"{gene_set}")
        print(f"Number of tissues with specific genes: {len(list(gene_set_dict_working_upset.index)[0])}")
        print(
            f"Number of genes specific to multiple tissues: {len([x for x in list(gene_set_dict_working_upset.index) if sum(x) > 1])}")

        gene_set_dict_working_upsetplot, gene_set_dict_working_upsetplot_axes = scaled_upset_plot(
            gene_set_dict_working_upset)
        gene_set_dict_working_upsetplot.savefig(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_genes_in_gene_sets_across_tissues_{gene_set}.svg",
            bbox_inches='tight')
        gene_set_dict_working_upsetplot.savefig(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_genes_in_gene_sets_across_tissues_{gene_set}.png",
            bbox_inches='tight')
        plt.close()

    return


def grab_gene_specificity_secretome_data(ensg_id: str, node_info: pd.DataFrame) -> Tuple[str | None, str | None, float]:
    specificity_data = node_info[node_info['node_id'] == ensg_id]['feature_specificity_category'].values[0]
    if ~isinstance(specificity_data, str):
        specificity_data = None
    secretome_data = node_info[node_info['node_id'] == ensg_id]['feature_secretome_location'].values[0]
    if ~isinstance(specificity_data, str):
        secretome_data = None
    expression_level = node_info[node_info['node_id'] == ensg_id]['feature_nTPM'].values[0]
    if not expression_level:
        expression_level = 0

    return specificity_data, secretome_data, expression_level


""" Gather metrics across all tissues """


def direct_overlaps(all_test_results: Dict[str, Any], analysis: str) -> pd.DataFrame:
    # read in olink proteins
    #  this mapping file comes from the construct_tissue_graphs pipeline
    olink_proteins = pd.read_csv(f"~/data/internal/proteomics/olink_field_uniprot_ensg_genesymbol_mapping.csv")
    high_missingness_proteins = pd.read_csv(f"~/data/internal/proteomics/high_missingness_protein_fields.txt")
    olink_post_qc_proteins = olink_proteins[
        ~olink_proteins.protein_field.isin(high_missingness_proteins.protein_field)].copy()
    olink_ensembl = sorted(list(olink_post_qc_proteins.ensg_id.unique()))

    all_direct_overlaps_list = []
    for tissue, tissue_dict in all_test_results.items():
        # read in relevant graph nodes info
        tissue_nodes = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/nodes.csv")
        for gene_set_name, gene_set_dict in tissue_dict.items():
            direct_overlap_result_from_analysis = gene_set_dict['direct_overlap']
            if direct_overlap_result_from_analysis['n_shared_genes'] == 0:
                continue
            direct_overlap_result = {'tissue': tissue, 'gene_set': gene_set_name,
                                     'n_shared': direct_overlap_result_from_analysis['n_shared_genes'],
                                     'pval': direct_overlap_result_from_analysis['fisher_greater_pval']}
            # how many of the shared genes (1) have some degree of tissue specificity (2) are secreted into the blood (3) both
            genes_with_specificity = set()
            genes_secreted_to_blood = set()
            genes_not_expressed_at_baseline = set()
            for gene in direct_overlap_result_from_analysis['shared_genes']:
                gene_specificity, gene_secretome, gene_ntpm = grab_gene_specificity_secretome_data(gene,
                                                                                                   node_info=tissue_nodes)
                if gene_ntpm < 0.1:
                    genes_not_expressed_at_baseline |= {gene}
                    continue
                if gene_specificity is not None:
                    genes_with_specificity |= {gene}
                if gene_secretome == 'Secreted to blood':
                    genes_secreted_to_blood |= {gene}
            genes_secreted_and_specific = genes_with_specificity & genes_secreted_to_blood - genes_not_expressed_at_baseline
            direct_overlap_result |= {
                'n_shared_at_baseline_exp': direct_overlap_result_from_analysis['n_shared_genes'] - len(
                    genes_not_expressed_at_baseline),
                'n_specific_genes': len(genes_with_specificity),
                'specific_genes': '-'.join(list(genes_with_specificity)),
                'n_blood_secreted_genes': len(genes_secreted_to_blood),
                'n_specific_secreted_exp': len(genes_secreted_and_specific)}
            # Reactome ORA for shared genes
            shared_genes_reactome_ora = reactome_ora(reactome_annots=reactome, tissue=tissue,
                                                     test_gene_set=direct_overlap_result_from_analysis['shared_genes'],
                                                     background=olink_ensembl)
            if not shared_genes_reactome_ora.empty:
                direct_overlap_result['shared_genes_reactome_ora'] = '|'.join(
                    shared_genes_reactome_ora['annotation'].to_list())
                direct_overlap_result['shared_genes_n_reactomes'] = len(shared_genes_reactome_ora)
            # Reactome ORA for tissue-specific shared genes
            if len(genes_with_specificity) > 0:
                shared_specific_genes_reactome_ora = reactome_ora(reactome_annots=reactome, tissue=tissue,
                                                                  test_gene_set=genes_with_specificity,
                                                                  background=olink_ensembl)
                if not shared_specific_genes_reactome_ora.empty:
                    direct_overlap_result['shared_specific_genes_reactome_ora'] = '|'.join(
                        shared_specific_genes_reactome_ora['annotation'].to_list())
                    direct_overlap_result['shared_genes_n_reactomes'] = len(shared_genes_reactome_ora)
            # store results
            all_direct_overlaps_list.append(direct_overlap_result)
    if len(all_direct_overlaps_list) == 0:
        return pd.DataFrame()
    all_direct_overlaps = pd.DataFrame(all_direct_overlaps_list)
    all_direct_overlaps = significance_labelling(all_direct_overlaps, pvalue_col='pval')

    return all_direct_overlaps


def get_shortest_paths(all_test_results: Dict[str, Any], analysis: str) -> pd.DataFrame:
    all_shortest_path_results = []
    for tissue, tissue_dict in all_test_results.items():
        for gene_set_name, gene_set_dict in tissue_dict.items():
            # cossim_shortest_path_from_analysis = gene_set_dict['shortest_paths_between_cossim_best_matches']
            shortest_path_from_analysis = gene_set_dict['shortest_paths_between_sets']
            all_shortest_path_results.append(pd.concat([
                # pd.DataFrame({'n_hops': ['-'.join(path_list).count('node') - 1 for path_list in
                #                          cossim_shortest_path_from_analysis.values()]})
                # .assign(tissue=tissue, gene_set=gene_set_name, path_type='Cossim'),
                pd.DataFrame({'n_hops': ['-'.join(path_list).count('node') - 1 for path_list in
                                         shortest_path_from_analysis.values()]})
                .assign(tissue=tissue, gene_set=gene_set_name, path_type='Topology')
            ]))
    if len(all_shortest_path_results) == 0:
        return pd.DataFrame()
    all_shortest_paths = pd.concat(all_shortest_path_results)

    return all_shortest_paths


""" Find shortest paths that have some tissue specificity """


def get_shortest_paths_with_specificity_secretion_info(all_test_results: Dict[str, Any], analysis: str):
    def annotation_categoriser(annotation_bool_vector: tuple[bool, bool, bool]) -> str:
        category_bool_vector_dict = {
            'Source only': (True, False, False),
            'Target only': (False, True, False),
            'Mediators only': (False, False, True),
            'Source & target': (True, True, False),
            'Source & mediators': (True, False, True),
            'Target & mediators': (False, True, True),
            'Source & target & mediators': (True, True, True)
        }
        bool_vector_category_dict = {v: k for k, v in category_bool_vector_dict.items()}

        return bool_vector_category_dict[annotation_bool_vector]

    shortest_paths_with_annotations_list = []
    shortest_paths_by_specificity_category_list = []  # tissue, gene_set, category, count, specificity-or-secretion
    shortest_paths_by_secretion_category_list = []
    for tissue, tissue_dict in all_test_results.items():
        tissue_nodes = pd.read_csv(f"~/data/internal/knowledge_graph/tissue_graphs/{tissue}/nodes.csv")
        tissue_specific_genes = tissue_nodes[tissue_nodes['feature_specificity_category'].notna()]['node_id'].to_list()
        secreted_genes = tissue_nodes[tissue_nodes['feature_secretome_location'] == 'Secreted to blood'][
            'node_id'].to_list()
        for gene_set_name, gene_set_dict in tissue_dict.items():
            specificity_category_dict = collections.defaultdict(int)
            secretion_category_dict = collections.defaultdict(int)
            for source_target, path in gene_set_dict['shortest_paths_between_sets'].items():
                node_list = [x.split(':')[1] for x in path if 'node' in x]
                if node_list[0] in tissue_specific_genes:
                    source_specific = 1
                else:
                    source_specific = 0
                if node_list[-1] in tissue_specific_genes:
                    target_specific = 1
                else:
                    target_specific = 0
                mediators_specific = list(set(node_list[1:-1]) & set(tissue_specific_genes))
                mediators_specific_count = len(set(node_list[1:-1]) & set(tissue_specific_genes))

                specificity_vector = (source_specific > 0, target_specific > 0, mediators_specific_count > 0)
                if specificity_vector != (False, False, False):
                    specificity_category = annotation_categoriser(specificity_vector)
                    specificity_category_dict[specificity_category] += 1

                if node_list[0] in secreted_genes:
                    source_secreted = 1
                else:
                    source_secreted = 0
                if node_list[-1] in secreted_genes:
                    target_secreted = 1
                else:
                    target_secreted = 0
                mediators_secreted = list(set(node_list[1:-1]) & set(secreted_genes))
                mediators_secreted_count = len(set(node_list[1:-1]) & set(secreted_genes))

                secretion_vector = (source_secreted > 0, target_secreted > 0, mediators_secreted_count > 0)
                if secretion_vector != (False, False, False):
                    secretion_category = annotation_categoriser(secretion_vector)
                    secretion_category_dict[secretion_category] += 1

                shortest_paths_with_annotations_list.append({
                    'tissue': tissue, 'gene_set': gene_set_name,
                    'source': source_target[0], 'target': source_target[1],
                    'length': len(node_list), 'mediators': '-'.join(node_list[1:-1]),
                    'source_specific_count': source_specific, 'target_specific_count': target_specific,
                    'mediators_specific_count': mediators_specific_count,
                    'mediators_specific': '-'.join(mediators_specific),
                    'source_secreted_count': source_secreted, 'target_secreted_count': target_secreted,
                    'mediators_secreted_count': mediators_secreted_count,
                    'mediators_secreted': '-'.join(mediators_secreted)
                })
            shortest_paths_by_specificity_category_list.append({'tissue': tissue, 'gene_set': gene_set_name} |
                                                               specificity_category_dict)
            shortest_paths_by_secretion_category_list.append({'tissue': tissue, 'gene_set': gene_set_name} |
                                                             secretion_category_dict)
    if len(shortest_paths_with_annotations_list) == 0:
        return pd.DataFrame(), pd.DataFrame()
    shortest_paths_with_annotations = pd.DataFrame(shortest_paths_with_annotations_list)
    shortest_paths_with_annotations_filtered = shortest_paths_with_annotations[
        shortest_paths_with_annotations['source_specific_count'] + shortest_paths_with_annotations[
            'target_specific_count'] + shortest_paths_with_annotations['mediators_specific_count'] > 0
        ].copy()

    # tally up types of path specificity/secretion info
    shortest_paths_by_specificity_category = pd.DataFrame(shortest_paths_by_specificity_category_list).fillna(0)
    shortest_paths_by_secretion_category = pd.DataFrame(shortest_paths_by_secretion_category_list).fillna(0)

    shortest_paths_by_annotation_category = pd.concat([
        shortest_paths_by_specificity_category.melt(id_vars=['tissue', 'gene_set'], value_name='counts',
                                                    var_name='category').assign(annotation_type='Specificity'),
        shortest_paths_by_secretion_category.melt(id_vars=['tissue', 'gene_set'], value_name='counts',
                                                  var_name='category').assign(annotation_type='Secretion')
    ])

    # plot
    # annotation_category_palette = {
    #     'All shared': 'tab:blue', 'Tissue-specific': 'tab:orange', 'Secreted to blood': 'tab:green',
    #                    'Secreted and specific': 'tab:red'}
    annotation_category_ordered = ['Source only', 'Target only', 'Mediators only', 'Source & target',
                                   'Source & mediators', 'Target & mediators', 'Source & target & mediators']
    shortest_paths_by_annotation_category_to_plot = shortest_paths_by_annotation_category.rename(
        columns={'gene_set': 'GWAS gene set',
                 'annotation_type': 'Annotation'})

    mpl.use('Agg')
    g = sns.FacetGrid(shortest_paths_by_annotation_category_to_plot, row='GWAS gene set', col='Annotation', height=8,
                      aspect=0.5, sharex=False)
    g.map_dataframe(sns.barplot, y='tissue', x='counts', hue='category', hue_order=annotation_category_ordered,
                    palette='tab10')
    g.set_axis_labels('Counts', 'Tissue')
    # g.tick_params(axis='x', rotation=90)
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    g.add_legend()
    for ax in g.axes_dict.values():
        ax.xaxis.grid(True, which='major', color='grey', alpha=0.6)
        ax.set_axisbelow(True)
    g.tight_layout()
    g.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/shortest_paths_with_specific_secretion_annotations.png")
    g.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/shortest_paths_with_specific_secretion_annotations.svg")
    plt.close()

    return shortest_paths_with_annotations_filtered, shortest_paths_by_annotation_category


def analyse_specific_shortest_paths(specificity_df: pd.DataFrame):
    # prioritise shortest paths that have specificity throughout their paths - then analyse, e.g. enrichment of all specific genes in the prioritised paths

    shortest_paths_with_specificity_secretion_info_ranked = specificity_df.copy()
    shortest_paths_with_specificity_secretion_info_ranked['specific_count_counter'] = (
                                                                                              shortest_paths_with_specificity_secretion_info_ranked[
                                                                                                  'source_specific_count'] > 0).astype(
        int) + (shortest_paths_with_specificity_secretion_info_ranked['target_specific_count'] > 0).astype(int) + (
                                                                                              shortest_paths_with_specificity_secretion_info_ranked[
                                                                                                  'mediators_specific_count'] > 0).astype(
        int)
    shortest_paths_with_max_specificity = shortest_paths_with_specificity_secretion_info_ranked[
        shortest_paths_with_specificity_secretion_info_ranked['specific_count_counter'] ==
        shortest_paths_with_specificity_secretion_info_ranked['specific_count_counter'].max()].sort_values(by='length')
    # account for if there are no mediators
    shortest_paths_with_max_specificity = pd.concat([shortest_paths_with_max_specificity,
                                                     shortest_paths_with_specificity_secretion_info_ranked[(
                                                                                                                   shortest_paths_with_specificity_secretion_info_ranked[
                                                                                                                       'length'] == 2) & (
                                                                                                                   shortest_paths_with_specificity_secretion_info_ranked[
                                                                                                                       'specific_count_counter'] == 2)]])
    print(
        f"Max shortest path 'specificity': {shortest_paths_with_specificity_secretion_info_ranked['specific_count_counter'].max()}")

    if shortest_paths_with_specificity_secretion_info_ranked['specific_count_counter'].max() < 2:
        print(f"Max shortest path 'specificity' is under 2: no paths are tissue-specific.")
        return pd.DataFrame()

    relevant_gene_sets = list(set(shortest_paths_with_max_specificity['gene_set']) - {'all_gwas'})
    for gene_set in relevant_gene_sets:
        # gene_set = 'pqtl'
        shortest_paths_with_max_specificity_specific_gene_set = shortest_paths_with_max_specificity[
            shortest_paths_with_max_specificity.gene_set == gene_set].fillna('').copy()

        # which tissues host paths of max specificity
        specific_shortest_paths_tissues = shortest_paths_with_max_specificity_specific_gene_set.tissue.value_counts()
        plt.figure(figsize=(8, 6))
        sns.barplot(specific_shortest_paths_tissues)
        plt.xlabel('Tissue')
        plt.ylabel('Number of max.-specificity paths')
        plt.tick_params(axis='x', rotation=90)
        plt.gca().yaxis.grid(True, which='major', color='grey', alpha=0.6)
        plt.gca().set_axisbelow(True)
        plt.tight_layout()
        plt.savefig(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_in_which_tissues_gwas-{gene_set}.svg")
        plt.savefig(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_in_which_tissues_gwas-{gene_set}.png")
        plt.close()

        genes_per_tissue = {
            tissue: set(tissue_df.source) | set(tissue_df.target) | set(
                '-'.join(list(tissue_df.mediators_specific)).split('-'))
            for tissue, tissue_df in shortest_paths_with_max_specificity_specific_gene_set.groupby('tissue')}

        # are path genes specific to each tissue unique to each tissue? (upset plot)
        if len(genes_per_tissue) > 1:
            genes_per_tissue_upset = upsetplot.from_contents({k: list(v) for k, v in genes_per_tissue.items()})
            genes_per_tissue_upsetplot, axes = scaled_upset_plot(genes_per_tissue_upset)
            genes_per_tissue_upsetplot.savefig(
                f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_gene_overlap_across_relevant_tissues_gwas-{gene_set}.svg",
                bbox_inches='tight')
            genes_per_tissue_upsetplot.savefig(
                f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_gene_overlap_across_relevant_tissues_gwas-{gene_set}.png",
                bbox_inches='tight')
            plt.close()

        # what are path genes specific to each tissue doing? (pathway enrichment)
        genes_per_tissue_reactome_ora = {}
        for tissue, tissue_gene_set in genes_per_tissue.items():
            genes_tissue_reactome_ora = reactome_ora(reactome, tissue=tissue, test_gene_set=tissue_gene_set)
            if genes_tissue_reactome_ora.empty:
                continue
            genes_tissue_reactome_ora = genes_tissue_reactome_ora.merge(
                reactome[['reactome_annotation_id', 'annotation_name']].drop_duplicates()
                .rename(columns={'reactome_annotation_id': 'annotation'}))
            genes_per_tissue_reactome_ora[tissue] = genes_tissue_reactome_ora
            genes_tissue_reactome_ora.to_csv(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/"
                                             f"specific_paths_gene_reactome_ora_{tissue}_gwas-{gene_set}.csv",
                                             index=False)
            _ = reactome_parentage_analysis(annotation_set=list(genes_tissue_reactome_ora.annotation.unique()),
                                            plot_path=f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/"
                                                      f"parentage_analysis_specific_paths_gene_reactome_ora_{tissue}_gwas-{gene_set}",
                                            export_plot=True)

    return shortest_paths_with_max_specificity


""" Execution """
ensembl_gene_set = generate_ensembl_gene_set(current_version='115')
ensg_to_symbol_map = dict(zip(ensembl_gene_set.gene_id, ensembl_gene_set.gene_name))

reactome_annotations_full = pd.read_csv(
    f"~/data/external/reactome/Ensembl2Reactome_PE_All_Levels.txt", sep='\t', header=None)
# prep annotation file
reactome = prep_annotation_files(raw_annotation_file=reactome_annotations_full, ontology='reactome')

mpl.use('Agg')
analyses = sorted([x.name for x in os.scandir('~/data/internal/knowledge_graph/gene_protein_graph_analysis') if x.is_dir()])
# test_analysis = 'gwasG30_coxG30'
for analysis in analyses:
    print(analysis)
    analysis_results_path = f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/graph_analysis_pickles"
    all_results = {}
    for result_file in [f"{analysis_results_path}/{x.name}" for x in os.scandir(analysis_results_path) if x.is_file()]:
        with open(result_file, 'rb') as pickle_file:
            all_results[result_file.split('/')[-1].split('.')[0]] = pickle.load(pickle_file)

    # make export dirs
    os.makedirs(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/", exist_ok=True)
    os.makedirs(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables/", exist_ok=True)

    # plot tissue specificity of GWAS and proteomics genes
    print('Results read. Plotting tissue specificity of the GWAS and proteomics gene sets.')
    plot_gene_set_specificity(all_test_results=all_results, analysis=analysis)

    # direct overlap
    print('Plotting complete. Conducting overlap analysis.')
    direct_overlap = direct_overlaps(all_test_results=all_results, analysis=analysis)
    if not direct_overlap.empty:
        direct_overlap.to_csv(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/direct_overlaps.csv",
                              index=False)

    # shortest path length
    print('Overlap analysis complete. Analysing and plotting shortest path lengths.')
    shortest_paths = get_shortest_paths(all_test_results=all_results, analysis=analysis)
    if not shortest_paths.empty:
        shortest_paths.to_csv(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/shortest_paths.csv",
                              index=False)

    print('Analysing shortest path lengths complete. Analysing shortest path specificity, and plotting.')
    # shortest path specificity
    (shortest_paths_with_specificity_secretion_info,
     shortest_paths_with_specificity_secretion_categories) = get_shortest_paths_with_specificity_secretion_info(
        all_test_results=all_results, analysis=analysis)
    if not shortest_paths_with_specificity_secretion_info.empty:
        shortest_paths_with_specificity_secretion_info.to_csv(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/"
            f"shortest_paths_with_specificity_secretion_info.csv",
            index=False)
        print('Prioritising paths based on specificity, and analysing.')
        specific_shortest_paths = analyse_specific_shortest_paths(
            specificity_df=shortest_paths_with_specificity_secretion_info)
        if not specific_shortest_paths.empty:
            specific_shortest_paths.to_csv(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/"
                                           f"{analysis}/tissue_specific_shortest_paths.csv",
                                           index=False)
            # make prioritised shortest paths tables with gene symbols for thesis
            specific_shortest_paths_symbols = specific_shortest_paths[
                specific_shortest_paths['gene_set'] != 'all_gwas'].sort_values(
                by=['tissue', 'gene_set'], ascending=True)
            specific_shortest_paths_symbols['source'] = specific_shortest_paths_symbols.source.map(ensg_to_symbol_map)
            specific_shortest_paths_symbols['target'] = specific_shortest_paths_symbols.target.map(ensg_to_symbol_map)
            specific_shortest_paths_symbols['mediators'] = specific_shortest_paths_symbols.mediators.apply(
                lambda x: '-'.join([ensg_to_symbol_map[gene] for gene in x.split('-')]) if pd.notna(x) else '')
            specific_shortest_paths_symbols[
                'mediators_specific'] = specific_shortest_paths_symbols.mediators_specific.apply(
                lambda x: '-'.join([ensg_to_symbol_map[gene] for gene in x.split('-')]) if pd.notna(x) else '')
            specific_shortest_paths_symbols['tissue'] = [
                f"{x.split('_')[0]} ({x.split('_')[1]} {' '.join([xx.lower() for xx in x.split('_')[2:]])})".replace(
                    ' )', ')')
                if '_' in x else x
                for x in specific_shortest_paths_symbols['tissue'].to_list()]
            specific_shortest_paths_symbols['tissue'] = specific_shortest_paths_symbols['tissue'].replace(
                'Small (Intestine terminal ileum)',
                'Small Intestine (Terminal ileum)')
            specific_shortest_paths_symbols[
                ['tissue', 'gene_set', 'source', 'source_specific_count', 'target', 'target_specific_count',
                 'mediators', 'mediators_specific']].to_csv(
                f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables/tissue_specific_shortest_paths_gene_symbols.csv",
                index=False)

    if not shortest_paths_with_specificity_secretion_categories.empty:
        shortest_paths_with_specificity_secretion_categories.to_csv(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/"
            f"shortest_paths_with_specificity_secretion_categories.csv",
            index=False)

    print(f"Analysis complete for {analysis}.")

# any overlap between shortest path enrichments between tissues? (upset)
for analysis in analyses:

    # read in csvs with 'specific_paths_gene_reactome_ora_*'
    reactome_results_dict = {
        reactome_file.name.replace('specific_paths_gene_reactome_ora_', '').replace('.csv', ''): pd.read_csv(
            f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/{reactome_file.name}")[
            'annotation'].to_list()
        for reactome_file in os.scandir(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/")  #
        if reactome_file.is_file() and 'specific_paths_gene_reactome_ora' in reactome_file.name}

    # set up upset
    if len(reactome_results_dict) == 1:
        continue
    reactomes_per_tissue_geneset_upset = upsetplot.from_contents(reactome_results_dict)
    fig = plt.figure()
    # totals_plot_elements=0 removes the horizontal per-set bar charts
    reactomes_per_tissue_geneset_upsetplot, axes = scaled_upset_plot(reactomes_per_tissue_geneset_upset)
    reactomes_per_tissue_geneset_upsetplot.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_reactome_overlap_across_tissues_genesets.svg",
        bbox_inches='tight')
    reactomes_per_tissue_geneset_upsetplot.savefig(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/specific_paths_reactome_overlap_across_tissues_genesets.png",
        bbox_inches='tight')
    plt.close()

    # export upset
    reactomes_per_tissue_geneset_upset.reset_index().merge(
        reactome[['reactome_annotation_id', 'annotation_name']]
        .rename(columns={'reactome_annotation_id': 'id'})).drop_duplicates().to_csv(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/"
        f"specific_paths_reactome_overlap_across_tissues_genesets.csv", index=False)

""" Select parentage analysis """
# T2D - OA
e11_oa_reactomes_upset = pd.read_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/gwasE11_coxM15-M16-M17-M18/"
    f"specific_paths_reactome_overlap_across_tissues_genesets.csv")
sets_of_interest = [
    ['Liver_gwas-magma', 'Liver_gwas-pqtl'],
    ['Spleen_gwas-magma', 'Spleen_gwas-pqtl'],
    ['Pancreas_gwas-pqtl', 'Liver_gwas-pqtl'],
    ['Breast_Mammary_Tissue_gwas-pqtl', 'Adipose_Visceral_Omentum_gwas-pqtl', 'Adipose_Subcutaneous_gwas-pqtl',
     'Adipose_Visceral_Omentum_gwas-magma']
]
for interest in sets_of_interest:
    interest_reactomes = set(e11_oa_reactomes_upset[e11_oa_reactomes_upset[interest[0]]]['id'])
    interest_reactomes = interest_reactomes.intersection(
        *[e11_oa_reactomes_upset[e11_oa_reactomes_upset[single_interest]]['id'].to_list()
          for single_interest in interest[1:]])
    _ = reactome_parentage_analysis(annotation_set=list(interest_reactomes),
                                    plot_path=f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/"
                                              f"gwasE11_coxM15-M16-M17-M18/plots/"
                                              f"parentage_analysis_specific_paths_reactome_overlap_parentage_{'_'.join(interest)}",
                                    export_plot=True)

# CAD - OA
i25_oa_reactomes_upset = pd.read_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/gwasI25_coxM15-M16-M17-M18/"
    f"specific_paths_reactome_overlap_across_tissues_genesets.csv")
sets_of_interest = [
    ['Liver_gwas-magma', 'Liver_gwas-pqtl'],
    ['Lung_gwas-pqtl', 'Spleen_gwas-pqtl'],
    ['Lung_gwas-pqtl', 'Colon_Sigmoid_gwas-pqtl'],
    ['Lung_gwas-pqtl', 'Testis_gwas-eqtl'],
    ['Testis_gwas-pqtl', 'Testis_gwas-eqtl'],
    ['Small_Intestine_Terminal_Ileum_gwas-pqtl', 'Spleen_gwas-pqtl']
]
for interest in sets_of_interest:
    interest_reactomes = set(i25_oa_reactomes_upset[i25_oa_reactomes_upset[interest[0]]]['id'])
    interest_reactomes = interest_reactomes.intersection(
        *[i25_oa_reactomes_upset[i25_oa_reactomes_upset[single_interest]]['id'].to_list()
          for single_interest in interest[1:]])
    _ = reactome_parentage_analysis(annotation_set=list(interest_reactomes),
                                    plot_path=f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/"
                                              f"gwasI25_coxM15-M16-M17-M18/plots/"
                                              f"parentage_analysis_specific_paths_reactome_overlap_parentage_{'_'.join(interest)}",
                                    export_plot=True)

# OA - cataracts
oa_cataracts_upset = pd.read_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/gwasM15-M16-M17-M18_coxH25-H26/"
    f"specific_paths_reactome_overlap_across_tissues_genesets.csv")
sets_of_interest = [
    ['Liver_gwas-eqtl', 'Liver_gwas-pqtl'],
    ['Liver_gwas-magma', 'Liver_gwas-eqtl', 'Liver_gwas-pqtl'],
    ['Spleen_gwas-magma', 'Spleen_gwas-pqtl'],
    ['Muscle_Skeletal_gwas-eqtl', 'Muscle_Skeletal_gwas-magma']
]
for interest in sets_of_interest:
    interest_reactomes = set(oa_cataracts_upset[oa_cataracts_upset[interest[0]]]['id'])
    interest_reactomes = interest_reactomes.intersection(
        *[oa_cataracts_upset[oa_cataracts_upset[single_interest]]['id'].to_list()
          for single_interest in interest[1:]])
    _ = reactome_parentage_analysis(annotation_set=list(interest_reactomes),
                                    plot_path=f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/"
                                              f"gwasM15-M16-M17-M18_coxH25-H26/plots/"
                                              f"parentage_analysis_specific_paths_reactome_overlap_parentage_{'_'.join(interest)}",
                                    export_plot=True)

""" Get pathways enriched in direct overlap genes """
reactome_annotations_full = pd.read_csv(
    f"graph_construction/database_data_files/reactome/Ensembl2Reactome_PE_All_Levels.txt", sep='\t', header=None)
# prep annotation file
reactome = prep_annotation_files(raw_annotation_file=reactome_annotations_full, ontology='reactome')
reactome_annotations_id_names = reactome[['reactome_annotation_id', 'annotation_name']].drop_duplicates(
    ignore_index=True)
reactome_id_to_name_map = dict(
    zip(reactome_annotations_id_names['reactome_annotation_id'], reactome_annotations_id_names['annotation_name']))

# THIS IS RUN FOR SPECIFIC ANALYSES
analyses = sorted([x.name for x in os.scandir('~/data/internal/knowledge_graph/gene_protein_graph_analysis') if x.is_dir()])
analysis = analyses[-2]
direct_overlap_results = pd.read_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/direct_overlaps.csv")
direct_overlap_results = direct_overlap_results[
    (direct_overlap_results['shared_genes_n_reactomes'] > 0) & (direct_overlap_results['fdr_sig'] == 1)].copy()

# make table for export
direct_overlap_table = direct_overlap_results[['tissue', 'gene_set', 'shared_genes_reactome_ora']].copy()
direct_overlap_table = direct_overlap_table[direct_overlap_table['gene_set'] != 'all_gwas'].copy()
direct_overlap_table['gene_set'] = direct_overlap_table['gene_set'].replace(
    {'eqtl': 'eQTL', 'pqtl': 'pQTL', 'magma': 'MAGMA'})
direct_overlap_table['tissue'] = [
    f"{x.split('_')[0]} ({x.split('_')[1]} {' '.join([xx.lower() for xx in x.split('_')[2:]])})".replace(' )', ')')
    if '_' in x else x
    for x in direct_overlap_table['tissue'].to_list()]
direct_overlap_table['tissue'] = direct_overlap_table['tissue'].replace('Small (Intestine terminal ileum)',
                                                                        'Small Intestine (Terminal ileum)')
direct_overlap_table['shared_genes_reactome_ora'] = ['; '.join([reactome_id_to_name_map[idd] for idd in x.split('|')])
                                                     for x in direct_overlap_table['shared_genes_reactome_ora']]
direct_overlap_table.to_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables/direct_overlap_reactome_ora_table.csv",
    index=False)

# plot raw numbers of annotations per tissue (facet by gene set)
direct_overlaps_reactome_counts = direct_overlap_results[['tissue', 'gene_set', 'shared_genes_n_reactomes']].copy()
direct_overlaps_reactome_counts = direct_overlaps_reactome_counts[
    direct_overlaps_reactome_counts['gene_set'] != 'all_gwas'].rename(columns={'gene_set': 'GWAS gene set'})
g = sns.FacetGrid(direct_overlaps_reactome_counts, col='GWAS gene set', height=7, aspect=0.4, sharex=False)
g.map_dataframe(sns.barplot, y='tissue', x='shared_genes_n_reactomes')
g.set_axis_labels('Counts', 'Tissue')
g.set_titles(col_template='{col_name}')
g.tick_params(axis='y', labelsize=8)
for ax in g.axes_dict.values():
    ax.xaxis.grid(True, which='major', color='grey', alpha=0.6)
    ax.set_axisbelow(True)
g.tight_layout()
g.savefig(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/n_reactomes_per_direct_overlap.png")
g.savefig(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/n_reactomes_per_direct_overlap.svg")
plt.close()

# plot upset of intersections (across all tissue-gene set combos)
direct_overlaps_reactome_sets = direct_overlap_results[['tissue', 'gene_set', 'shared_genes_reactome_ora']].copy()
direct_overlaps_reactome_sets = direct_overlaps_reactome_sets[direct_overlaps_reactome_sets['gene_set'] != 'all_gwas']
direct_overlaps_reactome_sets['tissue_gene_set'] = direct_overlaps_reactome_sets['tissue'] + '_gwas-' + \
                                                   direct_overlaps_reactome_sets['gene_set']

reactomes_per_overlap_upset = upsetplot.from_contents({k: v.split('|')
                                                       for k, v in
                                                       dict(zip(direct_overlaps_reactome_sets['tissue_gene_set'],
                                                                direct_overlaps_reactome_sets[
                                                                    'shared_genes_reactome_ora'])).items()})
reactomes_per_overlap_upsetplot, axes = scaled_upset_plot(reactomes_per_overlap_upset)
reactomes_per_overlap_upsetplot.savefig(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/direct_overlap_reactomes_per_tissue_gene_set.png")
reactomes_per_overlap_upsetplot.savefig(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/plots/direct_overlap_reactomes_per_tissue_gene_set.svg")
plt.close()

# export upset
reactomes_per_overlap_upset.reset_index().merge(
    reactome[['reactome_annotation_id', 'annotation_name']]
    .rename(columns={'reactome_annotation_id': 'id'})).drop_duplicates().to_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/"
    f"direct_overlap_reactome_overlap_across_tissues_genesets.csv", index=False)

""" Make thesis tables for Reactomes enriched in tissue-specific genes on tissue-specific paths """
analyses = sorted([x.name for x in os.scandir('~/data/internal/knowledge_graph/gene_protein_graph_analysis') if x.is_dir()])
analysis = analyses[5]
specific_path_enrichment_files = [x.name for x in
                                  os.scandir(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}")
                                  if x.is_file() and 'specific_paths_gene_reactome_ora' in x.name]
enrichment_df_dicts = []
for enrichment_file in specific_path_enrichment_files:
    enrichment_df = pd.read_csv(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/{enrichment_file}")
    info_name = enrichment_file.replace('specific_paths_gene_reactome_ora_', '')
    enrichment_tissue = info_name.split('_gwas')[0]
    enrichment_gene_set = info_name.split('_gwas-')[1].replace('.csv', '')
    enrichment_df_dicts.append({'tissue': enrichment_tissue, 'gene_set': enrichment_gene_set,
                                'reactomes': '; '.join(enrichment_df.annotation_name.to_list())})
enrichment_df_collated = pd.DataFrame(enrichment_df_dicts)
enrichment_df_collated['gene_set'] = enrichment_df_collated['gene_set'].replace(
    {'eqtl': 'eQTL', 'pqtl': 'pQTL', 'magma': 'MAGMA'})
enrichment_df_collated['tissue'] = [
    f"{x.split('_')[0]} ({x.split('_')[1]} {' '.join([xx.lower() for xx in x.split('_')[2:]])})".replace(' )', ')')
    if '_' in x else x
    for x in enrichment_df_collated['tissue'].to_list()]
enrichment_df_collated['tissue'] = enrichment_df_collated['tissue'].replace('Small (Intestine terminal ileum)',
                                                                            'Small Intestine (Terminal ileum)')
os.makedirs(f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables", exist_ok=True)
enrichment_df_collated.to_csv(
    f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables/specific_paths_gene_reactome_ora_table.csv",
    index=False)

""" Make thesis tables for prioritised paths with gene symbols """
ensembl_gene_set = generate_ensembl_gene_set('115')
ensg_to_symbol_map = dict(zip(ensembl_gene_set.gene_id, ensembl_gene_set.gene_name))

analyses = sorted([x.name for x in os.scandir('~/data/internal/knowledge_graph/gene_protein_graph_analysis') if x.is_dir()])
for analysis in analyses:
    prioritised_paths = pd.read_csv(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tissue_specific_shortest_paths.csv")
    prioritised_paths = prioritised_paths[prioritised_paths['gene_set'] != 'all_gwas'].sort_values(
        by=['tissue', 'gene_set'], ascending=True)
    prioritised_paths['source'] = prioritised_paths.source.map(ensg_to_symbol_map)
    prioritised_paths['target'] = prioritised_paths.target.map(ensg_to_symbol_map)
    prioritised_paths['mediators'] = prioritised_paths.mediators.apply(
        lambda x: '; '.join([ensg_to_symbol_map[gene] for gene in x.split('-')]) if pd.notna(x) else '')
    prioritised_paths['mediators_specific'] = prioritised_paths.mediators_specific.apply(
        lambda x: '; '.join([ensg_to_symbol_map[gene] for gene in x.split('-')]) if pd.notna(x) else '')
    prioritised_paths['tissue'] = [
        f"{x.split('_')[0]} ({x.split('_')[1]} {' '.join([xx.lower() for xx in x.split('_')[2:]])})".replace(' )', ')')
        if '_' in x else x
        for x in prioritised_paths['tissue'].to_list()]
    prioritised_paths['tissue'] = prioritised_paths['tissue'].replace('Small (Intestine terminal ileum)',
                                                                      'Small Intestine (Terminal ileum)')
    prioritised_paths['gene_set'] = prioritised_paths['gene_set'].replace({'pqtl': 'pQTL', 'eqtl': 'eQTL', 'magma': 'MAGMA'})
    prioritised_paths[
        ['tissue', 'gene_set', 'source', 'source_specific_count', 'target', 'target_specific_count', 'mediators',
         'mediators_specific']].to_csv(
        f"~/data/internal/knowledge_graph/gene_protein_graph_analysis/{analysis}/tables/tissue_specific_shortest_paths_gene_symbols.csv",
        index=False)
