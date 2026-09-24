from utils.significance_labelling import significance_labelling
from utils.enrichment_analyses import slim_annotations, propagation_up_ontology_hierarchy, overrepresentation

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.use('TkAgg')

""" Where do high-effect pathways from MAGMA localise to? """
# Reactome and GO annotation data
reactome_annotations_full = pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None)
reactome_annotations_full = reactome_annotations_full[
    (reactome_annotations_full[6] == 'TAS') & (reactome_annotations_full[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                     axis=1).set_axis(
    ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()
# propagation up hierarchy already dealt with by Reactome
# slimming not done in case slimming removes high-effect annotations found by MAGMA

# read in annotations file, name columns, filter, propagate up hierarchy, split by aspect, and slim
go_annotations = pd.read_csv('~/data/external/gene_ontology/goa_human.gaf.gz', sep='\t',
                          comment='!', header=None)
go_annotations.columns = ['id_db', 'id', 'symbol', 'relation', 'go_id', 'reference', 'evidence', 'additional_id',
                          'aspect', 'name', 'synonym', 'type', 'taxon', 'annotation_date', 'annotation_by',
                          'annotation_extension', 'gene_product_id']
go_annotations = go_annotations[(go_annotations['type'] == 'protein') &
                                (~go_annotations['evidence'].isin(['IEA', 'NAS', 'ND', 'NR'])) &
                                (~go_annotations['relation'].str.contains('NOT'))][
    ['id', 'go_id', 'aspect']].drop_duplicates(ignore_index=True)
# propagate annotations up hierarchy
go_annotations = propagation_up_ontology_hierarchy(id_ancestor_links=pd.read_csv(
    'graph_construction/database_data_files/gene_ontology/go_id_ancestor_relations.csv'),
    annotations_to_propagate=go_annotations,
    annotation_col_id='go_id')
# split by aspect and slim
go_annotations_by_aspect = {
    'go_' + aspect: slim_annotations(annotations=d, annotation_col_id='go_id', entity_col_id='id').rename(
        columns={'go_id': 'annotation_id'})
    for aspect, d in go_annotations.groupby('aspect')}

# get annotation names
annotation_names = pd.concat([reactome_annotations_full[['reactome_annotation_id', 'annotation_name']].rename(
    columns={'reactome_annotation_id': 'annotation'}),
                           pd.read_csv('~/data/external/gene_ontology/go_id_names_aspects.csv')[
                               ['id', 'name']].rename(
                               columns={'id': 'annotation', 'name': 'annotation_name'})]).drop_duplicates(
    ignore_index=True)

# MAGMA
magma_enrichment_results = pd.read_csv('~/data/internal/genomics/magma/lifestyles/reactome_enrichment_results_full.csv')
magma_enrichment_high_effect = {disease: d.sort_values('BETA', ascending=False, ignore_index=True).loc[:19,
                                         ['VARIABLE', 'BETA']].set_index('VARIABLE')['BETA'].to_dict()
                                for disease, d in magma_enrichment_results.groupby('disease')}
high_effect_genes_all_diseases = set()
for disease, disease_dict in magma_enrichment_high_effect.items():
    # identify Olink proteins annotated with the top 20 high-effect annotations for each disease
    high_effect_annots_genes = reactome_annotations_full[
        reactome_annotations_full['reactome_annotation_id'].isin(list(disease_dict.keys()))]['id'].to_list()
    # collect proteins
    high_effect_genes_all_diseases = high_effect_genes_all_diseases | set(high_effect_annots_genes)

high_effect_annots_genes_all_diseases_cc_enrichment = overrepresentation(test_set=list(high_effect_genes_all_diseases),
                                                                         background_annotations=go_annotations_by_aspect['go_C'],
                                                                         entity_col_id='id',
                                                                         annotation_col_id='annotation_id')
high_effect_annots_genes_all_diseases_cc_enrichment_labelled = significance_labelling(
    high_effect_annots_genes_all_diseases_cc_enrichment)
high_effect_annots_genes_all_diseases_cc_enrichment_labelled = high_effect_annots_genes_all_diseases_cc_enrichment_labelled.pd.merge(annotation_names)
high_effect_annots_genes_all_diseases_cc_enrichment_fdr_sig = \
high_effect_annots_genes_all_diseases_cc_enrichment_labelled[
    high_effect_annots_genes_all_diseases_cc_enrichment_labelled['fdr_sig'] == 1].pd.merge(annotation_names)
high_effect_annots_genes_all_diseases_cc_enrichment_labelled.to_csv('~/data/internal/genomics/magma/lifestyles/high_effect_pathways_gocc_enrichment.csv', index=False)

# FDR-sig GO-CC
go_c_fdr_sig_list = high_effect_annots_genes_all_diseases_cc_enrichment_fdr_sig['annotation_name'].to_list()

# plot
magma_cc_to_plot = pd.read_csv('~/data/internal/genomics/lifestyles/high_effect_pathways_gocc_enrichment.csv')
magma_cc_to_plot_sig = magma_cc_to_plot[magma_cc_to_plot['fdr_sig'] == 1].copy()
theme_assignment = {
    # Membrane receptor & channel complexes
    'Bcl-2 family protein complex': 'Membrane receptor &\n channel complexes',
    'receptor complex': 'Membrane receptor &\n channel complexes',
    'plasma membrane signaling receptor complex': 'Membrane receptor &\n channel complexes',
    'voltage-gated potassium channel complex': 'Membrane receptor &\n channel complexes',
    'potassium channel complex': 'Membrane receptor &\n channel complexes',
    'acetylcholine-gated channel complex': 'Membrane receptor &\n channel complexes',
    # Endocytic / vesicular trafficking
    'endocytic vesicle membrane': 'Endocytic /\nvesicular trafficking',
    'coated vesicle membrane': 'Endocytic /\nvesicular trafficking',
    'clathrin-coated vesicle membrane': 'Endocytic /\nvesicular trafficking',
    'clathrin-coated endocytic vesicle membrane': 'Endocytic /\nvesicular trafficking',
    'clathrin-coated endocytic vesicle': 'Endocytic /\nvesicular trafficking',
    'endolysosome': 'Endocytic /\nvesicular trafficking',
    'endolysosome membrane': 'Endocytic /\nvesicular trafficking',
    'lysosomal lumen': 'Endocytic /\nvesicular trafficking',
    # Nuclear / chromatin
    'telomere cap complex': 'Nuclear /\nchromatin',
    'nuclear telomere cap complex': 'Nuclear /\nchromatin',
    'protein-DNA complex': 'Nuclear /\nchromatin',
    'nucleosome': 'Nuclear /\nchromatin',
    # Extracellular / circulating
    'blood microparticle': 'Extracellular /\ncirculating',
    'immunoglobulin complex': 'Extracellular /\ncirculating',
    'immunoglobulin complex, circulating': 'Extracellular /\ncirculating',
}

magma_cc_to_plot_sig['Theme'] = magma_cc_to_plot_sig['annotation_name'].map(theme_assignment)
magma_cc_to_plot_sig['Theme'] = magma_cc_to_plot_sig.Theme.fillna('Misc.')
magma_cc_to_plot_sig['gene_ratio'] = magma_cc_to_plot_sig.n_test_annotated / magma_cc_to_plot_sig.n_total_annotated
magma_cc_to_plot_sig = magma_cc_to_plot_sig.sort_values(by='gene_ratio', ascending=False)
magma_cc_to_plot_sig['annotation_name'] = pd.Categorical(magma_cc_to_plot_sig['annotation_name'], magma_cc_to_plot_sig['annotation_name'].to_list())

fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=magma_cc_to_plot_sig, x='gene_ratio', y='annotation_name', hue='Theme', palette='tab10',
            hue_order=['Extracellular /\ncirculating', 'Membrane receptor &\n channel complexes', 'Endocytic /\nvesicular trafficking', 'Nuclear /\nchromatin', 'Misc.'])
ax.set_ylabel('Gene Ontology (Cellular Component)')
ax.set_xlabel('Gene ratio')
ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
plt.tight_layout()
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/high_effect_magma_enriched_go_cc_terms.svg')
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/high_effect_magma_enriched_go_cc_terms.png')
plt.close()
