import pandas as pd
import matplotlib as mpl
import numpy as np
import math

from utils.enrichment_analyses import prep_annotation_files
from utils.semantic_similarity import plot_high_flying_micas


mpl.use('TkAgg')

association = 'Firth'  # 'Cox'

""" Identify top-contributing MICAs """


def count_high_flying_micas(selected_micas_df: pd.DataFrame) -> pd.DataFrame:
    high_flying_micas_slim = selected_micas_df[['mica', 'disease1', 'disease2']].drop_duplicates()
    high_flying_micas_diseases = [
        {'mica': mica,
         'diseases': '|'.join([field_to_code_chapter_mapper[x]
                               for x in
                               sorted(list(set(mica_diseases_df['disease1']) | set(mica_diseases_df['disease2'])))])}
        for mica, mica_diseases_df in high_flying_micas_slim.groupby('mica')]
    high_flying_micas_diseases_df = pd.DataFrame(high_flying_micas_diseases)
    high_flying_micas_diseases_df['n_diseases'] = high_flying_micas_diseases_df['diseases'].apply(
        lambda x: len(x.split('|')))
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        pd.DataFrame(high_flying_micas_slim['mica'].value_counts().rename('pair_count')).reset_index())

    # get reactome pathway names
    reactome_pathways = pd.read_csv('graph_construction/database_data_files/reactome/ReactomePathways.txt', sep='\t',
                                    header=None, names=['id', 'name', 'species'])
    reactome_pathways = reactome_pathways[reactome_pathways['species'] == 'Homo sapiens'].drop(columns='species')
    high_flying_micas_diseases_df = high_flying_micas_diseases_df.merge(
        reactome_pathways.rename(columns={'id': 'mica', 'name': 'mica_name'}))

    return high_flying_micas_diseases_df.sort_values(by='pair_count', ascending=False).reset_index(drop=True)


# aux info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
field_to_code_chapter_mapper = dict(zip(disease_info.disease_field, disease_info.code_chapter))
bma_mica_contributions = pd.read_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_disease_resnik_micas.csv")

# protein background
olink_uniprot_mappings = pd.read_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv')
high_missingness_proteins = pd.read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/proteomics/high_missingness_protein_fields.txt')
all_qc_proteins = olink_uniprot_mappings[
    ~olink_uniprot_mappings['protein_field'].isin(high_missingness_proteins['protein_field'])].copy()
protein_background = all_qc_proteins['uniprot'].to_list()

# annotation info
annotations = prep_annotation_files(pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None),
    ontology='reactome').rename(columns={'reactome_annotation_id': 'annotation_id'})[
    ['id', 'annotation_id']].drop_duplicates(ignore_index=True)
hierarchy = pd.read_csv(
    '~/data/external/reactome/reactome_id_ancestor_relations.csv').set_axis(
    ['id', 'ancestors'], axis=1)
annotations_for_protein_background = annotations[annotations['id'].isin(protein_background)].copy()
annotations_for_protein_background = annotations_for_protein_background[['id', 'annotation_id']].drop_duplicates()

annotations_ic = -np.log(annotations_for_protein_background.annotation_id.value_counts() / len(
    annotations_for_protein_background.id.unique()))

""" Driver pathways """

# top 5 per disease pair by contribution + MICA IC >= 5.587705 (25th percentile of all Reactome ICs)
top5_micas_per_disease_pair = pd.concat(
    [disease_pair_df.iloc[:5, :] for _, disease_pair_df in bma_mica_contributions.groupby(['disease1', 'disease2'])])
prioritised_micas = top5_micas_per_disease_pair[
    top5_micas_per_disease_pair['resnik'] >= annotations_ic.quantile(0.25)].copy()
prioritised_micas_counted = count_high_flying_micas(prioritised_micas)

# look at these MICAs but index by how many disease-disease connections they're involved in given the set of diseases they're linked to
#  so if a MICA is linked to 10 diseases, there's a maximum of 45 pairs - how many of these is the MICA representing/linked to
prioritised_micas_counted['pairwise_completion'] = [
    prioritised_micas_counted.at[i, 'pair_count'] / math.comb(prioritised_micas_counted.at[i, 'n_diseases'],
                                                              2) for i in
    range(len(prioritised_micas_counted))]
(prioritised_micas_counted.pairwise_completion == 1).sum()
prioritised_micas_complete_connectors = prioritised_micas_counted[
    prioritised_micas_counted['pairwise_completion'] == 1].drop(columns='pairwise_completion').copy()
prioritised_micas_complete_connectors = prioritised_micas_complete_connectors[
    prioritised_micas_complete_connectors['n_diseases'] > 2].copy()
prioritised_micas_counted.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_semsim_common_driver_pathways.csv",
    index=False)
prioritised_micas_complete_connectors.to_csv(
    f"~/data/internal/proteomics/pathway_based_similarities/{association}_semsim_common_driver_pathways_complete_connectors.csv",
    index=False)

# plot
plot_path = f"~/ch3_proteomics/3.6 pathway_based_similarities/plots/{association}_semsim_common_driver_pathways"
plot_high_flying_micas(counted_micas_df=prioritised_micas_complete_connectors, n_top_pathways=25,
                       plot_path=plot_path)

# parentage analysis
from utils.enrichment_analyses import reactome_parentage_analysis

_ = reactome_parentage_analysis(annotation_set=prioritised_micas_complete_connectors.mica.to_list(),
                                plot_path=f"ch3_proteomics/3.6 pathway_based_similarities/plots/{association}_semsim_common_driver_pathways_parentage",
                                export_plot=True)
