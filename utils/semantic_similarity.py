""" Functions for pathway-based semantic similarity work. """
from typing import Tuple, List
import numpy as np
import pandas as pd
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


# calculate similarity metric (resnik, lin, jiang-conrath) for a given set of nodes
def weighted_semantic_similarity(node_id_weights_dict, ic_data, hierarchy_data, metric='resnik'):
    # node_id_weights_dict = {<node>: <weight>, ...}
    # ic_data = pd.Series with index = node_id, values = information content
    # hierarchy_data = pd.DataFrame with child node_id and all ancestor node_id for that child node_id in long format
    common_ancestors = set.intersection(
        *[set(hierarchy_data[hierarchy_data['id'] == x]['ancestors']) for x in node_id_weights_dict.keys()])
    if not bool(common_ancestors):
        return 0

    mean_weights = sum(list(node_id_weights_dict.values())) / len(node_id_weights_dict)

    mica = list(ic_data[list(common_ancestors)].sort_values(ascending=False).index)[0]
    resnik = ic_data[mica]

    if metric == 'resnik':
        similarity = resnik
    elif metric == 'lin':
        similarity = resnik * (len(node_id_weights_dict) / sum(ic_data[list(node_id_weights_dict.keys())]))
    elif metric == 'jiang-conrath':  # this is a distance metric, so invert to make a similarity metric
        similarity = 1 / sum(ic_data[list(node_id_weights_dict.keys())]) - (len(node_id_weights_dict) * resnik)
    else:
        print(f"{metric} is not covered. Returning None")
        return None

    return similarity * mean_weights


def semantic_similarity(node_ids, ic_data, hierarchy_data, metric='resnik'):
    """
    Calculate similarity metric (resnik, lin, jiang-conrath) for a given set of nodes
    :param node_ids:
    :param ic_data:
    :param hierarchy_data:
    :param metric:
    :return:
    """
    # node_ids = [<node>, <node>, ...]
    # ic_data = pd.Series with index = node_id, values = information content
    # hierarchy_data = pd.DataFrame with child node_id and all ancestor node_id for that child node_id in long format
    common_ancestors = set.intersection(
        *[set(hierarchy_data[hierarchy_data['id'] == x]['ancestors']) for x in node_ids])
    if not bool(common_ancestors):
        return 0, 'none'

    mica = list(ic_data[list(common_ancestors)].sort_values(ascending=False).index)[0]
    resnik = ic_data[mica]

    if metric == 'resnik':
        similarity = resnik
    elif metric == 'lin':
        similarity = resnik * (len(node_ids) / sum(ic_data[list(node_ids)]))
    elif metric == 'jiang-conrath':  # this is a distance metric, so invert to make a similarity metric
        similarity = 1 / sum(ic_data[list(node_ids)]) - (len(node_ids) * resnik)
    else:
        print(f"{metric} is not covered. Returning None")
        return None, None

    return similarity, mica


def calculate_mica_contributions(sem_sim_array: np.array, mica_array: np.array) -> pd.DataFrame:
    """
    Calculate contribution of individual MICAs to a Resnik BMA score between two annotation sets.
    :param sem_sim_array:
    :param mica_array:
    :return:
    """
    # one-liner: identify MICAs driving a Resnik BMA score between two annotation sets
    # given storage matrices for semantic similarity and corresponding MICAs
    # for each annotation in setA (row indices), identify best-match MICA via max(Resnik)
    # then the same for setB (col indices)
    # set max(Resnik)s into list + rank by contribution to the overall Resnik BMA score

    # identify which MICA is responsible for each element in each set's max(Resnik) (i.e., best match)
    set_a_maxes = sem_sim_array.max(axis=1)  # max Resnik per row
    set_a_maxes_idxs = sem_sim_array.argmax(axis=1)  # col idx of max per row
    set_a_max_micas = mica_array[np.arange(mica_array.shape[0]), set_a_maxes_idxs]  # MICA with max Resnik per row

    set_b_maxes = sem_sim_array.max(axis=0)  # max Resnik per col
    set_b_maxes_idxs = sem_sim_array.argmax(axis=0)  # row idx of max per col
    set_b_max_micas = mica_array[set_b_maxes_idxs, np.arange(mica_array.shape[1])]  # MICA with max Resnik per col

    # quantify 'contribution': (BMA_full - BMA_loo) / BMA_full. BMA_loo = BMA with one MICA (and its Resnik) removed
    # rank by this
    baseline_bma = 0.5 * (set_a_maxes.mean() + set_b_maxes.mean())
    set_a_loo_bma = 0.5 * ((set_a_maxes.sum() - set_a_maxes) / (set_a_maxes.size - 1) + set_b_maxes.mean())
    set_a_maxes_contributions = (baseline_bma - set_a_loo_bma) / baseline_bma
    set_b_loo_bma = 0.5 * ((set_b_maxes.sum() - set_b_maxes) / (set_b_maxes.size - 1) + set_a_maxes.mean())
    set_b_maxes_contributions = (baseline_bma - set_b_loo_bma) / baseline_bma

    mica_df = pd.concat([pd.DataFrame({'resnik': set_a_maxes, 'mica': set_a_max_micas,
                                       'contribution': set_a_maxes_contributions}),
                         pd.DataFrame({'resnik': set_b_maxes, 'mica': set_b_max_micas,
                                       'contribution': set_b_maxes_contributions})])

    return mica_df.sort_values(by='contribution', ascending=False).drop_duplicates()


def resnik_bma_between_two_annotation_sets(
        annot1: List[str], annot2: List[str], annotation_info_content: pd.DataFrame, annotation_hierarchy: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate Resnik BMA, including a measure of individual MICA contribution, between two annotation sets
    :param annot1:
    :param annot2:
    :param annotation_info_content:
    :param annotation_hierarchy:
    :return:
    """
    # disease_annotation_list_of_dicts = [disease_annotation_weights[disease] for disease in disease_group]

    # set up sem-sim storage array
    sem_sim_storage_array = np.zeros(shape=(len(annot1), len(annot2)))
    # sem_sim_storage_array = np.zeros(shape=tuple([len(x) for x in disease_annotation_list_of_dicts]))
    # set up MICA storage array
    mica_storage_array = np.empty(shape=(len(annot1), len(annot2)), dtype=object)

    # set up annotation-to-arraycol mappings
    disease_annotation_maps_to_arraycols = [{e: i for i, e in enumerate(annot1)}, {e: i for i, e in enumerate(annot2)}]
    # disease_annotation_maps_to_arraycols = [dict(zip(list(x.keys()), range(len(x)))) for x in
    #                                         disease_annotation_list_of_dicts]

    annotation_combinations = list(itertools.product(annot1, annot2))
    # annotation_combinations = list(itertools.product(*[list(x.keys()) for x in disease_annotation_list_of_dicts]))
    for annotation_combination in annotation_combinations:
        annotation_combination_similarity = semantic_similarity(
            node_ids=annotation_combination,
            ic_data=annotation_info_content,
            hierarchy_data=annotation_hierarchy)
        # store result in storage array
        sem_sim_storage_array[tuple(disease_annotation_maps_to_arraycols[i][x] for i, x in
                                    enumerate(annotation_combination))] = annotation_combination_similarity[0]
        mica_storage_array[tuple(disease_annotation_maps_to_arraycols[i][x] for i, x in
                                 enumerate(annotation_combination))] = annotation_combination_similarity[1]

    # perform best-match averaging - max(axis=1) covers setA, max(axis=0) covers setB
    bma_sem_sim = 0.5 * (sem_sim_storage_array.max(axis=1).mean() + sem_sim_storage_array.max(axis=0).mean())

    # semsim result df
    sem_sim_result = pd.DataFrame([{'bma_sem_sim': bma_sem_sim,
                                    'n_annots_disease1': len(annot1),
                                    'n_annots_disease2': len(annot2)}])

    # get MICA contributions and store result
    mica_contributions = calculate_mica_contributions(sem_sim_array=sem_sim_storage_array, mica_array=mica_storage_array)

    return sem_sim_result, mica_contributions


def get_min_max_resnik_bma_for_reactome():
    """ Max, min SemSim for a given ontology (Reactome in this case) """
    from utils.enrichment_analyses import prep_annotation_files
    import numpy as np
    annotations = prep_annotation_files(pd.read_csv(
        'graph_construction/database_data_files/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None),
        ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
        ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
    # hierarchy = pd.read_csv(
    #     'graph_construction/database_data_files/reactome/reactome_id_ancestor_relations.csv').set_axis(
    #     ['id', 'ancestors'], axis=1)
    annotations_ic = -np.log(annotations.annotation_id.value_counts() / len(annotations.id.unique()))
    upper_info_content_bound = -np.log(10 / 19026)
    # there will be no terms used in SemSim calcs that are too specific (because they will have been removed in the enrichment analysis step), so remove them
    # but there will be terms that are too general to be used in enrichment used in SemSim calcs, so retain these
    # hence filter only on upper_info_content_bound
    annotations_ic_no_specifics = annotations_ic[annotations_ic < upper_info_content_bound]
    # max SemSim is if all terms tested are identical and are the lowest-level terms possible (mean of the 20 highest ICs)
    resnik_bma_max = annotations_ic_no_specifics.sort_values(ascending=False).values[:20].mean()

    # min SemSim is if all terms tested exist so far apart there is no MICA for any combination of them => 0
    print(f"Max. Resnik BMA: {resnik_bma_max} \nMin. Resnik BMA: 0")


def plot_high_flying_micas(counted_micas_df: pd.DataFrame, n_top_pathways: int, plot_path: str):
    """
    Given a set of high-flying driver MICAs, plot their links with diseases and disease pairs via a
    (1) heatmap (MICA-disease)
    (2) bar chart counting n(diseases)
    (3) bar chart counting n(disease pairs)
    The MICAs are ordered by n(disease pairs) because we are interested in MICAs that drive similarity between pairs.
    :param counted_micas_df:
    :param n_top_pathways:
    :param plot_name:
    :return:
    """

    disease_info = pd.read_csv(
        'ukbiobank/ard_identification/by_age_of_onset_stats/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
        str) + ')'

    # heatmap. mica_name on y-axis (cut-off after n characters), ARDs x-axis.
    # Separate bar charts of n_diseases and pair_count, as with the combined manhattan-ish plot
    counted_micas_df = counted_micas_df.reset_index(drop=True)

    # split mica_name over 2 lines if n_characters > 50
    def mica_name_split(name):
        if len(name) <= 50:
            return name

        lines = []
        current_line = ""
        words = name.split(" ")

        for word in words:
            # If adding this word would exceed the limit
            if current_line and len(current_line) + 1 + len(word) > 50:
                lines.append(current_line)
                current_line = word
            else:
                if current_line:
                    current_line += " " + word
                else:
                    current_line = word

        if current_line:
            lines.append(current_line)

        return "\n".join(lines)

    counted_micas_df['mica_name'] = counted_micas_df['mica_name'].map(mica_name_split)

    # convert counted_micas_df to heatmap format
    counted_micas_heatmap = pd.DataFrame(index=counted_micas_df['mica_name'], columns=disease_info['code_chapter'],
                                         dtype=int)

    for _, row in counted_micas_df.iterrows():
        counted_micas_heatmap.loc[row['mica_name'], row['diseases'].split('|')] = 1
    counted_micas_heatmap.fillna(0, inplace=True)

    top_n_pathways_counts = counted_micas_df.iloc[:n_top_pathways, :].copy()

    # Figure layout
    fig = plt.figure(figsize=(16, 12))
    fig.set_constrained_layout(True)
    gs = fig.add_gridspec(1, 3, width_ratios=[8, 1, 1], wspace=0.05)

    # --- Scatter plot ---
    ax = fig.add_subplot(gs[0])

    sns.heatmap(
        data=counted_micas_heatmap.iloc[:n_top_pathways, :], linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
        cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
        linecolor='black'
    )

    ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
    ax.set_ylabel('Reactome term')
    ax.set_xticklabels(ax.get_xticklabels(), size=10)
    ax.set_yticklabels(ax.get_yticklabels(), size=12)

    # --- Bar chart (counts per disease) ---
    # log_vars = np.log1p(n_vars_per_disease_full.values)
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    ax2.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['n_diseases'], color='grey',
             alpha=0.7)
    ax2.set_xlabel('n(diseases)')
    ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

    # Add vertical grid lines at tick marks
    xticks = ax2.get_xticks()
    for xt in xticks:
        ax2.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['n_diseases'].to_list()):
        ax2.text(count, y, str(count), va='center', fontsize=8)

    # --- Bar chart (gene counts per disease) ---
    # log_genes = np.log1p(n_genes_per_disease_full.values)
    ax3 = fig.add_subplot(gs[2], sharey=ax)
    ax3.barh(y=[x + 0.5 for x in top_n_pathways_counts.index], width=top_n_pathways_counts['pair_count'], color='grey',
             alpha=0.7)
    ax3.set_xlabel('n(disease pairs)')
    ax3.tick_params(axis='y', left=False, labelleft=False)

    # Vertical grid lines
    for xt in ax3.get_xticks():
        ax3.axvline(x=xt, color='lightgrey', lw=0.8, ls='-', zorder=0)

    # Add raw counts labels
    for y, count in zip([x + 0.5 for x in top_n_pathways_counts.index], top_n_pathways_counts['pair_count'].to_list()):
        ax3.text(count, y, str(count), va='center', fontsize=8)

    plt.show()

    plt.savefig(f"{plot_path}.png")
    plt.savefig(f"{plot_path}.svg")

    plt.close()

    return


