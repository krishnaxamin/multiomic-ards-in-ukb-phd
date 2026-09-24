""" Functions to conduct annotation enrichment analyses. """
import numpy as np
import typing
import pandas as pd
import scipy


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


def propagation_up_ontology_hierarchy(id_ancestor_links, annotations_to_propagate, annotation_col_id,
                                      relation_class=None, relation_class_col_id=None):
    # argument check
    if relation_class is not None:
        assert relation_class_col_id is not None, 'To specify a relation_class, the col_id must also be specified.'
    if relation_class is not None:
        id_ancestor_links = id_ancestor_links[id_ancestor_links[relation_class_col_id] == relation_class].copy()
    else:
        id_ancestor_links = id_ancestor_links.copy()

    annotations_to_propagate_list = list(annotations_to_propagate[annotation_col_id].unique())
    filtered_links = id_ancestor_links[id_ancestor_links['id'].isin(annotations_to_propagate_list)].copy()

    # generate mapping: {id: [list of ancestors, including id itself]}
    id_to_ancestor_map = (filtered_links
                          .groupby('id')['ancestor']
                          .unique()
                          .apply(list)
                          .to_dict())

    annotations_to_propagate_with_ancestors = annotations_to_propagate.copy()
    annotations_to_propagate_with_ancestors['ancestors'] = annotations_to_propagate_with_ancestors[
        annotation_col_id].apply(
        lambda x: id_to_ancestor_map[x])

    annotations_propagated = (annotations_to_propagate_with_ancestors.explode('ancestors', ignore_index=True)
                              .drop(columns=annotation_col_id)
                              .rename(columns={'ancestors': annotation_col_id})
                              .drop_duplicates(ignore_index=True)
                              .merge(annotations_to_propagate.assign(original_annotation=1),
                                     on=list(annotations_to_propagate.columns),
                                     how='outer')
                              .fillna(0)
                              .astype({'original_annotation': 'int8'}))

    return annotations_propagated


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


def prep_annotation_files(raw_annotation_file, ontology):
    if ontology == 'reactome':
        prepped_annotation_file = raw_annotation_file[
            (raw_annotation_file[6] == 'TAS') & (raw_annotation_file[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                 axis=1).set_axis(
            ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()
    elif ontology == 'geneontology':
        raw_annotation_file.columns = ['id_db', 'id', 'symbol', 'relation', 'go_id', 'reference', 'evidence',
                                       'additional_id',
                                       'aspect', 'name', 'synonym', 'type', 'taxon', 'annotation_date', 'annotation_by',
                                       'annotation_extension', 'gene_product_id']
        go_annotations = raw_annotation_file[(raw_annotation_file['type'] == 'protein') &
                                             (~raw_annotation_file['evidence'].isin(['IEA', 'NAS', 'ND', 'NR'])) &
                                             (~raw_annotation_file['relation'].str.contains('NOT'))][
            ['id', 'go_id', 'aspect']].drop_duplicates(ignore_index=True)
        # propagate annotations up hierarchy
        prepped_annotation_file = propagation_up_ontology_hierarchy(id_ancestor_links=pd.read_csv(
            'graph_construction/database_data_files/gene_ontology/go_id_ancestor_relations.csv'),
            annotations_to_propagate=go_annotations,
            annotation_col_id='go_id')
    else:
        print(f"{ontology} is not covered.")
        return

    return prepped_annotation_file


def reactome_parentage_analysis(annotation_set: typing.List[str], plot_path: str = None, export_plot: bool = False):
    """
    Given a set of Reactome annotations, identify all ancestors and quantify the proportion of terms that are
    descendants of each ancestor ('popularity'). Test whether 'popularity' is greater than expected, then export and
    plot these data.
    :param annotation_set:
    :param plot_path:
    :param export_plot:
    :return:
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    from utils.significance_labelling import significance_labelling

    mpl.use('TkAgg')

    id_ancestor_links = pd.read_csv(
        'graph_construction/database_data_files/reactome/reactome_id_ancestor_relations.csv')
    id_ancestor_links = id_ancestor_links[id_ancestor_links.id != id_ancestor_links.ancestors].copy()

    # get number of descendants covered by each ancestor of the annotation set
    # popularity = proportion of annotation set covered by an ancestor
    ancestor_popularity = pd.DataFrame(
        id_ancestor_links[id_ancestor_links.id.isin(annotation_set)].value_counts('ancestors')).reset_index().set_axis(
        ['annotation', 'n_descendants'], axis=1)
    ancestor_popularity['popularity'] = ancestor_popularity.n_descendants / len(annotation_set)

    # get reactome pathway names
    reactome_pathways = pd.read_csv('~/data/external/reactome/ReactomePathways.txt', sep='\t',
                                    header=None, names=['id', 'name', 'species'])
    reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
    ancestor_popularity = ancestor_popularity.merge(reactome_pathways.set_axis(['annotation', 'reactome'], axis=1))

    # what proportion of all Reactome terms are descendant of a given term
    ancestor_popularity_baseline = pd.DataFrame(
        id_ancestor_links[id_ancestor_links.ancestors.isin(ancestor_popularity.annotation)].value_counts(
            'ancestors') / id_ancestor_links.id.nunique()).reset_index().set_axis(['annotation', 'baseline_popularity'],
                                                                                  axis=1)
    ancestor_popularity = ancestor_popularity.merge(ancestor_popularity_baseline)

    # binomial test to find significantly greater popularity than expected
    ancestor_popularity_pval = []
    for row in ancestor_popularity.to_dict(orient='records'):
        ancestor_popularity_pval.append(scipy.stats.binomtest(k=row['n_descendants'],
                                                              n=len(annotation_set),
                                                              p=row['baseline_popularity'],
                                                              alternative='greater').pvalue)
    ancestor_popularity['pval'] = ancestor_popularity_pval

    # label with FDR < 0.05
    ancestor_popularity = significance_labelling(ancestor_popularity)
    ancestor_popularity = ancestor_popularity[['annotation', 'popularity', 'reactome', 'fdr_sig']].copy()

    # set up for plot
    ancestor_popularity_to_plot = ancestor_popularity.copy()
    ancestor_popularity_to_plot['popularity'] = ancestor_popularity_to_plot.popularity * 100
    ancestor_popularity_to_plot.loc[ancestor_popularity_to_plot.fdr_sig == 1, 'fdr_sig'] = 'Yes'
    ancestor_popularity_to_plot.loc[ancestor_popularity_to_plot.fdr_sig == 0, 'fdr_sig'] = 'No'

    # dynamic sizing of plot
    BASE_N = 27
    BASE_FIG_W = 8
    BASE_FIG_H = 6
    # BASE_FONT  = 12
    BASE_TICK = 10
    BASE_LABEL = 12

    scale = len(ancestor_popularity_to_plot) / BASE_N

    fig_w = max(8, BASE_FIG_W * (scale ** 0.25))  # width grows slower than height since terms are listed on y-axis
    fig_h = max(6, BASE_FIG_H * (scale ** 0.5))
    # font = max(8, BASE_FONT * (scale ** -0.3))  # fonts shrink as n grows
    # fonts shrink as n grows, but if n shrinks, cap font size at default so font doesn't get too large
    if scale <= 1:
        tick = BASE_TICK
        label = BASE_LABEL
    else:
        tick = max(5, BASE_TICK * (scale ** -0.3))  # fonts shrink as n grows, but if n shrinks, cap font size
        label = max(6, BASE_LABEL * (scale ** -0.2))  # axis label shrinks more slowly

    plt.figure(figsize=(fig_w, fig_h))
    sns.barplot(data=ancestor_popularity_to_plot, x='popularity', y='reactome', hue='fdr_sig', hue_order=['Yes', 'No'],
                palette=[plt.cm.coolwarm(1.0), 'lightgrey'])
    plt.xlabel('% terms for which X is an ancestor', fontsize=label)
    if scale < 1:
        plt.xlabel('% terms for which X \nis an ancestor', fontsize=label)
    plt.ylabel('Ancestor', fontsize=label)
    plt.legend(title='FDR < 0.05')
    plt.gca().tick_params(axis='both', labelsize=tick)
    plt.gca().xaxis.grid(True, which='major', color='grey', alpha=0.2)
    plt.gca().set_axisbelow(True)
    plt.tight_layout()

    if export_plot and plot_path is not None:
        plt.savefig(f"{plot_path}.svg")
        plt.savefig(f"{plot_path}.png")
        plt.close()

    return ancestor_popularity


def effect_vs_logp_plot(enrichment_results_df: pd.DataFrame, effect_size_col: str, plot_export_path: str,
                        pval_col: str = 'pval'):
    """
    Produce effect size-vs-log(p) plot to examine whether enrichment analyses are well-powered. Poorly powered analyses
    show high effect-high pval density.
    :param enrichment_results_df:
    :param effect_size_col:
    :param plot_export_path:
    :param pval_col:
    :return:
    """

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns

    mpl.use('Agg')

    enrichment_working = enrichment_results_df.copy()

    enrichment_working['logp'] = -1 * np.log10(enrichment_working[pval_col])
    enrichment_working.loc[enrichment_working[effect_size_col] == np.inf, effect_size_col] = round(enrichment_working[effect_size_col].max() / 100) * 100
    enrichment_working['log_or'] = np.log10(abs(enrichment_working[effect_size_col]))

    plt.figure(figsize=(8, 6))
    sns.scatterplot(enrichment_working, x='logp', y=effect_size_col)
    plt.xlabel('-log10(p)')
    plt.ylabel('Effect size')
    plt.axhline(y=0, color='grey', linestyle='dashed')
    plt.tight_layout()
    plt.savefig(f"{plot_export_path}.png")
    plt.savefig(f"{plot_export_path}.svg")
    plt.close()


def enrichment_qqplots(enrichment_results_df: pd.DataFrame, disease_col: str, disease_col_disease_id_type: str,
                       plot_export_path: str, pval_col: str = 'pval'):
    """
    Plot QQ-plots for each ARD. Little-to-no or negative deviation from the diagonal indicates underpowered analyses.
    :param enrichment_results_df:
    :param disease_col:
    :param disease_col_disease_id_type:
    :param plot_export_path:
    :param pval_col:
    :return:
    """

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns

    mpl.use('Agg')

    enrichment_working = enrichment_results_df.copy()

    disease_info = pd.read_csv('ukbiobank/ard_identification/by_age_of_onset_stats/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
        str) + ')'

    if disease_col_disease_id_type != 'code_chapter':
        enrichment_working[disease_col] = enrichment_working[disease_col].map(dict(zip(disease_info[disease_col_disease_id_type],
                                                                                       disease_info['code_chapter'])))

    def prep_qq_plot_data(pvalues):
        """ Given a set of pvalues, obtain the observed and expected log10p arrays for QQ plots"""
        # log10 transform - these are observed log10 values
        log10_pvalues = -1 * np.log10(pvalues)
        # arrange descending (highest log10p at top)
        log10_pvalues_sorted = np.sort(log10_pvalues)[::-1]
        # add column ('rank') going from 1 to n (highest log10p has '1', lowest has n)
        rank = np.array(range(1, (len(log10_pvalues_sorted) + 1)))
        # expected log10 values = -1 * np.log10(rank / (len(pvalues) + 1))
        expected_log10 = -1 * np.log10(rank / (len(pvalues) + 1))
        qqplot_data = pd.DataFrame({'obs_log10': log10_pvalues_sorted, 'exp_log10': expected_log10})
        return qqplot_data

    # set up figure for QQ plots
    figure, axis = plt.subplots(ncols=7, nrows=10, figsize=[11.7, 16.5])  # A3 size
    axis = axis.flatten()

    for disease_idx, disease in enumerate(list(disease_info['code_chapter'])):

        disease_enrichment = enrichment_working[enrichment_working[disease_col] == disease].copy()
        disease_enrichment_qqplot_data = prep_qq_plot_data(disease_enrichment[pval_col])

        sns.scatterplot(data=disease_enrichment_qqplot_data, x='exp_log10', y='obs_log10', edgecolor='none', alpha=0.6,
                        ax=axis[disease_idx])
        axis[disease_idx].axline((0, 0), slope=1, color='grey', linestyle='dashed')
        axis[disease_idx].set_title(disease)
        axis[disease_idx].set(xticklabels=[], yticklabels=[], ylabel=None, xlabel=None)
        axis[disease_idx].tick_params(bottom=False, left=False)

    # remove boxes 69 and 70 from the arrayed plot
    for i in range(1, 3):
        axis[7 * 10 - i].spines['top'].set_visible(False)
        axis[7 * 10 - i].spines['right'].set_visible(False)
        axis[7 * 10 - i].spines['bottom'].set_visible(False)
        axis[7 * 10 - i].spines['left'].set_visible(False)

        sns.histplot([], ax=axis[7 * 10 - i])
        axis[7 * 10 - i].set(xlabel=None, ylabel=None)
        axis[7 * 10 - i].tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    figure.tight_layout()
    figure.savefig(f"{plot_export_path}.png")
    figure.savefig(f"{plot_export_path}.svg")
    plt.close()
