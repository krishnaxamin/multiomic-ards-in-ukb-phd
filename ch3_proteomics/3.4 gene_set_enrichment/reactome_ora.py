""" Bottom-up Reactome enrichment analyses of FDR-sig associated proteins, without first-degree interactors. """
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyprind import ProgBar

from utils.significance_labelling import significance_labelling
from utils.enrichment_analyses import overrepresentation, slim_annotations

import numpy as np
import sys

mpl.use('TkAgg')

""" Define task """
association = 'firth'

""" Read invariant data """
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Read in ontology data """
reactome_annotations_full = pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None)
reactome_annotations_full = reactome_annotations_full[
    (reactome_annotations_full[6] == 'TAS') & (reactome_annotations_full[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                     axis=1).set_axis(
    ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()
# propagation up hierarchy already dealt with by Reactome

""" Read in protein data """
proteomics = pd.read_csv(f"~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_{association}.csv")
sig_proteomics = proteomics[proteomics['fdr_sig'] == 1].copy()

""" Read in mapping data """
# read in protein -> UniProt conversions
olink_uniprot_mappings = pd.read_csv('~/data/internal/proteomics/olink_field_uniprot_genesymbol_mapping.csv')
sig_proteomics_uniprot = sig_proteomics.merge(
    olink_uniprot_mappings[['uniprot', 'protein_field']].rename(columns={'protein_field': 'protein'}))

""" Define protein and annotation background """
# --- protein background ---
# use the mapping file, as this contains all the UKB proteins before QC
# remove the 3 proteins calculated as high missingness during QC
high_missingness_proteins = pd.read_csv('~/ch1_ard_discovery/1.2 cohort_data/1.2.1 qc/proteomics/high_missingness_protein_fields.txt')
all_qc_proteins = olink_uniprot_mappings[
    ~olink_uniprot_mappings['protein_field'].isin(high_missingness_proteins['protein_field'])].copy()
protein_background = all_qc_proteins['uniprot'].to_list()

# --- annotation background ---
annotations_for_protein_background = reactome_annotations_full[
    reactome_annotations_full['id'].isin(protein_background)].copy()
annotations_for_protein_background = annotations_for_protein_background[
    ['id', 'reactome_annotation_id', 'annotation_name']].drop_duplicates()

# slim
annotations_for_protein_background = slim_annotations(annotations=annotations_for_protein_background,
                                                      annotation_col_id='reactome_annotation_id', entity_col_id='id')

""" Annotate all assoc proteins """
sig_proteomics_with_annotations = annotations_for_protein_background[
    annotations_for_protein_background['id'].isin(sig_proteomics_uniprot.uniprot.to_list())].copy()

# rudimentary check that slimming has worked
assert (-np.log(sig_proteomics_with_annotations.reactome_annotation_id.value_counts() / len(
    sig_proteomics_with_annotations.id.unique()))).max() < -np.log(10 / 19026)

""" Overrepresentation enrichment """
overrep_results = pd.DataFrame()
bar = ProgBar(len(sig_proteomics.disease.unique()), stream=sys.stdout, title='Overrepresentation enrichment')
for disease in sig_proteomics.disease.unique():
    protein_test_set = sig_proteomics_uniprot[sig_proteomics_uniprot.disease == disease].uniprot.to_list()

    overrep_result = overrepresentation(test_set=protein_test_set,
                                        background_annotations=annotations_for_protein_background,
                                        entity_col_id='id', annotation_col_id='reactome_annotation_id')
    overrep_results = pd.concat([overrep_results, overrep_result.assign(disease=disease)])

    bar.update()

# tidy
overrep_results = (overrep_results.reset_index(drop=True)
                   .merge(reactome_annotations_full[['reactome_annotation_id', 'annotation_name']]
                          .rename(columns={'reactome_annotation_id': 'annotation'}))
                   .merge(disease_info[['disease_field', 'icd10_three_letter']]
                          .rename(columns={'disease_field': 'disease'}))).drop_duplicates(ignore_index=True)

# significance labelling
overrep_results_labelled = significance_labelling(overrep_results, fdr_permissive_group_by='disease', alpha=0.05)
overrep_results_nom_sig = overrep_results_labelled[overrep_results_labelled['nom_sig'] == 1].copy()

# export
overrep_results.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora.csv",
    index=False)
overrep_results_nom_sig.to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora_nom_sig.csv",
    index=False)

""" Plot """

# read in data
overrep_results_nom_sig = pd.read_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora_nom_sig.csv")
overrep_results_sig = overrep_results_nom_sig[overrep_results_nom_sig.fdr_sig == 1].copy()
overrep_results_sig = overrep_results_sig[['annotation', 'disease', 'annotation_name']].copy()

disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
overrep_results_sig = overrep_results_sig.merge(
    disease_info[['disease_field', 'code_chapter']].rename(columns={'disease_field': 'disease'}))

# pathways per disease
pathways_per_disease = pd.DataFrame(overrep_results_sig.value_counts('code_chapter')).reset_index()
plt.figure(figsize=(8, 6))
sns.barplot(data=pathways_per_disease, x='code_chapter', y='count')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('Number of FDR < 0.05 Reactome annotations')
plt.gca().tick_params(axis='x', rotation=90)
for container in plt.gca().containers:
    plt.gca().bar_label(container)
plt.gca().grid(True, which='major', axis='y', color='lightgrey', lw=0.8, ls='--')
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig(
    f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/n_reactomes_per_disease_{association}.png")
plt.savefig(
    f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/n_reactomes_per_disease_{association}.svg")
plt.close()

# diseases per pathway
diseases_per_pathway = overrep_results_sig.value_counts('annotation_name')

# pathways that are associated with multiple diseases
multidisease_pathways = overrep_results_sig[
    overrep_results_sig.annotation_name.isin(list(diseases_per_pathway[diseases_per_pathway > 1].index))][
    ['code_chapter', 'annotation_name']]
multidisease_pathways_heatmap = multidisease_pathways.assign(mark=1).pivot(index='annotation_name',
                                                                           columns='code_chapter',
                                                                           values='mark').fillna(0)

if association == 'firth':
    n_disease_limit = 5
    x_tick_fontsize = 10
    y_tick_fontsize = 10
elif association == 'cox':
    n_disease_limit = 45
    x_tick_fontsize = 8
    y_tick_fontsize = 9
else:
    raise ValueError('Wrong algorithm')

# Figure layout
fig = plt.figure(figsize=(10, 8))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)

# --- Scatter plot ---
ax = fig.add_subplot(gs[0])

sns.heatmap(
    data=multidisease_pathways_heatmap, linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
    cbar=False, cmap=mpl.colors.LinearSegmentedColormap.from_list('white_to_red', ['white', plt.cm.coolwarm(1.0)]),
    linecolor='black'
)

ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Reactome')
ax.set_xticklabels(ax.get_xticklabels(), size=7, rotation=90)
ax.set_yticklabels(ax.get_yticklabels(), size=9)

# --- Bar chart (counts per disease) ---
# log_vars = np.log1p(n_vars_per_disease_full.values)
multidisease_pathways_ax2 = pd.DataFrame(diseases_per_pathway[diseases_per_pathway > 1]).reset_index().set_axis(
    ['reactome', 'freq'], axis=1).sort_values(by='reactome').reset_index(drop=True)
ax2 = fig.add_subplot(gs[1], sharey=ax)
ax2.barh(y=[x + 0.5 for x in multidisease_pathways_ax2.index], width=multidisease_pathways_ax2['freq'], color='grey',
         alpha=0.7)
ax2.set_xlabel('n(diseases)')
ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

# Add vertical grid lines at tick marks
# ax2.set_xticks([2, 4, 6, 8])
ax2.grid(True, axis='x', color='lightgrey', lw=0.8, ls='--')
ax2.set_axisbelow(True)

plt.savefig(f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/diseases_per_reactome_with_many_diseases.png")
plt.savefig(f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/diseases_per_reactome_with_many_diseases.svg")
plt.close()

""" Parentage analysis """
from utils.enrichment_analyses import reactome_parentage_analysis

_ = reactome_parentage_analysis(annotation_set=list(overrep_results_sig.annotation.unique()),
                                plot_path=f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/parentage_analysis_all_sig_terms_{association}",
                                export_plot=True)

""" Disease-pathway table """
results_collated_list = []
for _, df in overrep_results_sig.groupby('code_chapter'):
    results_collated_list.append(df.assign(reactomes='; '.join(df.annotation_name.to_list())).drop(
        columns=['annotation', 'annotation_name', 'disease']).drop_duplicates())
results_collated = pd.concat(results_collated_list)
results_collated = results_collated.merge(disease_info[['icd10_three_letter', 'code_chapter']])
results_collated[['icd10_three_letter', 'reactomes']].to_csv(
    f"~/data/internal/proteomics/enrichment_analysis/{association}/disease_reactome_mappings.csv", index=False)

""" Power analyses """
from utils.enrichment_analyses import effect_vs_logp_plot, enrichment_qqplots
enrichment_results = pd.read_csv(f"~/data/internal/proteomics/enrichment_analysis/{association}/reactome_ora.csv")

# effect-vs-logp plot
effect_vs_logp_plot(enrichment_results_df=enrichment_results,
                    effect_size_col='odds_ratio', pval_col='pval',
                    plot_export_path=f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/reactome_enrichment_effect_vs_logp_plot")

# QQ plots
enrichment_qqplots(enrichment_results_df=enrichment_results,
                   disease_col='disease', disease_col_disease_id_type='disease_field', pval_col='pval',
                   plot_export_path=f"~/ch3_proteomics/3.4 gene_set_enrichment/{association}/plots/reactome_enrichment_qqplots")
