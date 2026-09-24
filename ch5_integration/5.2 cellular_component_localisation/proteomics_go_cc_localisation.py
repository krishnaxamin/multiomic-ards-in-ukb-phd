""" Test the Olink proteins for enrichment relative to all annotated proteins. """
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from utils.significance_labelling import significance_labelling
from utils.enrichment_analyses import overrepresentation, slim_annotations, propagation_up_ontology_hierarchy

mpl.use('TkAgg')

""" Read in data """
all_proteins = pd.read_csv('~/data/external/proteomics/olink_assay.dat', sep='\t')

# remove the 3 proteins calculated as high missingness during QC
high_missingness_proteins = pd.read_csv('~/data/internal/proteomics/high_missingness_protein_fields.txt')
all_qc_proteins = all_proteins[
    ~all_proteins['Assay'].isin(high_missingness_proteins['protein_field'].str.upper())].copy()

all_proteins_expanded = all_qc_proteins.assign(UniProt=all_qc_proteins['UniProt'].str.split('_')).explode('UniProt',
                                                                                                          ignore_index=True)
all_proteins_list = all_proteins_expanded['UniProt'].to_list()

annotations_dict = {}

""" Read in Reactomes """
reactome_annotations_full = pd.read_csv(
    '~/data/external/reactome/UniProt2Reactome_PE_All_Levels.txt', sep='\t', header=None)
reactome_annotations_full = reactome_annotations_full[
    (reactome_annotations_full[6] == 'TAS') & (reactome_annotations_full[7] == 'Homo sapiens')].drop([4, 6, 7],
                                                                                                     axis=1).set_axis(
    ['id', 'reactome_id', 'reactome_name', 'reactome_annotation_id', 'annotation_name'], axis=1).copy()
# propagation up hierarchy already dealt with by Reactome
# slim full background annotation
reactome_annotations_slimmed = slim_annotations(annotations=reactome_annotations_full,
                                                annotation_col_id='reactome_annotation_id', entity_col_id='id')

annotations_dict['reactome'] = reactome_annotations_slimmed.rename(columns={'reactome_annotation_id': 'annotation_id'})

""" Read in Gene Ontology """
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
    '~/data/external/gene_ontology/go_id_ancestor_relations.csv'),
    annotations_to_propagate=go_annotations,
    annotation_col_id='go_id')
# split by aspect and slim
go_annotations_by_aspect = {
    'go_' + aspect: slim_annotations(annotations=d, annotation_col_id='go_id', entity_col_id='id').rename(
        columns={'go_id': 'annotation_id'})
    for aspect, d in go_annotations.groupby('aspect')}

annotations_dict.update(go_annotations_by_aspect)

""" Assemble annotation names """
annotation_names = pd.concat([reactome_annotations_full[['reactome_annotation_id', 'annotation_name']].rename(
    columns={'reactome_annotation_id': 'annotation'}),
    pd.read_csv('~/data/external/gene_ontology/go_id_names_aspects.csv')[
        ['id', 'name']].rename(
        columns={'id': 'annotation', 'name': 'annotation_name'})]).drop_duplicates(
    ignore_index=True)

""" Overrepresentation enrichment """
overrep_results = pd.DataFrame()
for annot, annot_df in annotations_dict.items():
    overrep_result = overrepresentation(test_set=all_proteins_list, background_annotations=annot_df,
                                        entity_col_id='id', annotation_col_id='annotation_id')
    overrep_result_labelled = significance_labelling(overrep_result)
    overrep_results = pd.concat([overrep_results, overrep_result_labelled.assign(ontology=annot)])

# add names
overrep_results_named = overrep_results.merge(annotation_names)
overrep_results_nom_sig = overrep_results_named[overrep_results_named['pval'] < 0.05].copy()

# export
overrep_results_named.to_csv('~/data/internal/proteomics/enrichment_analysis/olink_vs_all.csv', index=False)
overrep_results_nom_sig.to_csv('~/data/internal/proteomics/enrichment_analysis/olink_vs_all_nom_sig.csv', index=False)

# count number of FDR-sig enrichments per ontology/aspect
overrep_results[overrep_results['fdr_sig'] == 1].ontology.value_counts()

""" Is their enrichment in cellular compartments indicating bias towards secretion? """
overrep_results_named = pd.read_csv('~/data/internal/proteomics/enrichment_analysis/olink_vs_all.csv')
# FDR-sig GO-CC
go_c_fdr_sig_list = \
    overrep_results_named[(overrep_results_named['ontology'] == 'go_C') & (overrep_results_named['fdr_sig'] == 1)][
        'annotation_name'].to_list()
go_c_fdr_sig_df = overrep_results_named[
    (overrep_results_named['ontology'] == 'go_C') & (overrep_results_named['fdr_sig'] == 1)].copy()
go_c_fdr_sig_olink = go_annotations[(go_annotations['go_id'].isin(go_c_fdr_sig_df['annotation'])) & (
    go_annotations['id'].isin(all_proteins_list))].copy()
go_c_fdr_sig_olink_unique_proteins = go_c_fdr_sig_olink.id.unique().tolist()
# the 31 FDR-sig GO_CC pathways are broadly localised in Secretory / blood circulatory annotations, e.g. Secretory vesiclbes, plasma HDL

# plot these terms with an assignment of whether they are: Intracellular; Secretory; Extracellular
go_c_fdr_sig_df = overrep_results_named[
    (overrep_results_named['ontology'] == 'go_C') & (overrep_results_named['fdr_sig'] == 1)].copy()
theme_assignment = {
    # Secretory
    'specific granule': 'Secretory',
    'tertiary granule': 'Secretory',
    'secretory granule': 'Secretory',
    'azurophil granule': 'Secretory',
    'platelet alpha granule': 'Secretory',
    'platelet dense granule': 'Secretory',
    'Golgi lumen': 'Secretory',
    'platelet alpha granule lumen': 'Secretory',
    'vesicle lumen': 'Secretory',
    'Secretory granule lumen': 'Secretory',
    'cytoplasmic vesicle lumen': 'Secretory',
    'endoplasmic reticulum lumen': 'Secretory',
    'specific granule lumen': 'Secretory',
    'tertiary granule lumen': 'Secretory',
    'azurophil granule lumen': 'Secretory',
    'platelet dense granule lumen': 'Secretory',
    # Extracellular
    'external encapsulating structure': 'Extracellular /\ncirculating',
    'extracellular matrix': 'Extracellular /\ncirculating',
    'external side of plasma membrane': 'Extracellular /\ncirculating',
    'plasma lipoprotein particle': 'Extracellular /\ncirculating',
    'high-density lipoprotein particle': 'Extracellular /\ncirculating',
    'lipoprotein particle': 'Extracellular /\ncirculating',
    'protein-lipid complex': 'Extracellular /\ncirculating',
    'immunological synapse': 'Extracellular /\ncirculating',
    'protein complex involved in cell adhesion': 'Extracellular /\ncirculating',
    # Intracellular
    'receptor complex': 'Intracellular',
    'A band': 'Intracellular',
    'vacuolar lumen': 'Intracellular - degradation',
    'lysosomal lumen': 'Intracellular - degradation',
}

go_c_fdr_sig_df['Theme'] = go_c_fdr_sig_df['annotation_name'].map(theme_assignment)
go_c_fdr_sig_df['Theme'] = go_c_fdr_sig_df.Theme.fillna('Misc.')
go_c_fdr_sig_df['gene_ratio'] = go_c_fdr_sig_df.n_test_annotated / go_c_fdr_sig_df.n_total_annotated
go_c_fdr_sig_df = go_c_fdr_sig_df.sort_values(by='gene_ratio', ascending=False)
go_c_fdr_sig_df['annotation_name'] = pd.Categorical(go_c_fdr_sig_df['annotation_name'],
                                                    go_c_fdr_sig_df['annotation_name'].to_list())

fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=go_c_fdr_sig_df, x='gene_ratio', y='annotation_name', hue='Theme', palette='tab20',
            hue_order=['Extracellular /\ncirculating', 'Secretory', 'Intracellular', 'Misc.'])
ax.set_ylabel('Gene Ontology (Cellular Component)')
ax.set_xlabel('Gene ratio')
ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
plt.tight_layout()
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/olink_vs_all_go_cc_terms.svg')
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/olink_vs_all_go_cc_terms.png')
plt.close()

""" Do high-effect pathways localise to Secretory / circulatory annotations? """
disease_enrichment_results = pd.read_csv(
    '~/data/internal/proteomics/enrichment_analysis/cox/reactome_ora.csv')
disease_enrichment_high_effect = {disease: d.sort_values('odds_ratio', ascending=False, ignore_index=True).loc[:19,
                                           ['annotation', 'odds_ratio']].set_index('annotation')['odds_ratio'].to_dict()
                                  for disease, d in disease_enrichment_results.groupby('disease')}
# for each disease, find Olink proteins annotated with the top 20 annotations, and test for GO-CC enrichment
# high_effect_annots_olinks_cc_enrichment_all_diseases = pd.DataFrame()
high_effect_annots_olinks_all_diseases = set()
for disease, disease_dict in disease_enrichment_high_effect.items():
    # identify Olink proteins annotated with the top 20 high-effect annotations for each disease
    high_effect_annots_olinks = reactome_annotations_full[(reactome_annotations_full['id'].isin(all_proteins_list)) & (
        reactome_annotations_full['reactome_annotation_id'].isin(list(disease_dict.keys())))]['id'].to_list()
    # collect proteins
    high_effect_annots_olinks_all_diseases = high_effect_annots_olinks_all_diseases | set(high_effect_annots_olinks)

high_effect_annots_olinks_all_diseases_cc_enrichment = overrepresentation(
    test_set=list(high_effect_annots_olinks_all_diseases),
    background_annotations=go_annotations_by_aspect['go_C'],
    entity_col_id='id',
    annotation_col_id='annotation_id')
high_effect_annots_olinks_all_diseases_cc_enrichment_labelled = significance_labelling(
    high_effect_annots_olinks_all_diseases_cc_enrichment)
high_effect_annots_olinks_all_diseases_cc_enrichment_named = high_effect_annots_olinks_all_diseases_cc_enrichment_labelled.merge(
    annotation_names)
high_effect_annots_olinks_all_diseases_cc_enrichment_fdr_sig = \
    high_effect_annots_olinks_all_diseases_cc_enrichment_labelled[
        high_effect_annots_olinks_all_diseases_cc_enrichment_labelled['fdr_sig'] == 1].merge(annotation_names)
high_effect_annots_olinks_all_diseases_cc_enrichment_named.to_csv(
    '~/data/internal/proteomics/enrichment_analysis/high_effect_cox_pathways_gocc_enrichment.csv', index=False)

# plot
high_effect_cc_to_plot = pd.read_csv(
    '~/data/internal/proteomics/enrichment_analysis/high_effect_cox_pathways_gocc_enrichment.csv')
high_effect_cc_to_plot_sig = high_effect_cc_to_plot[high_effect_cc_to_plot['fdr_sig'] == 1].copy()
theme_assignment = {
    # Secretory
    'Golgi lumen': 'Secretory',
    'vesicle lumen': 'Secretory',
    'tertiary granule': 'Secretory',
    'secretory granule': 'Secretory',
    'secretory granule lumen': 'Secretory',
    'cytoplasmic vesicle lumen': 'Secretory',
    'endoplasmic reticulum lumen': 'Secretory',
    'specific granule': 'Secretory',
    'specific granule lumen': 'Secretory',
    'coated vesicle membrane': 'Secretory',
    'tertiary granule lumen': 'Secretory',
    'ER to Golgi transport vesicle membrane': 'Secretory',
    'COPII-coated ER to Golgi transport vesicle': 'Secretory',
    'multivesicular body': 'Secretory',
    'platelet alpha granule': 'Secretory',
    'platelet alpha granule lumen': 'Secretory',
    'azurophil granule lumen': 'Secretory',
    'endoplasmic reticulum-Golgi intermediate compartment': 'Secretory',
    'endoplasmic reticulum-Golgi intermediate compartment membrane': 'Secretory',
    'transport vesicle membrane': 'Secretory',
    # Extracellular / circulating
    'external side of plasma membrane': 'Extracellular / circulating',
    'extrinsic component of plasma membrane': 'Extracellular / circulating',
    'basement membrane': 'Extracellular / circulating',
    'external encapsulating structure': 'Extracellular / circulating',
    'extracellular matrix': 'Extracellular / circulating',
    'specialized extracellular matrix': 'Extracellular / circulating',
    'protein-lipid complex': 'Extracellular / circulating',
    'plasma lipoprotein particle': 'Extracellular / circulating',
    'very-low-density lipoprotein particle': 'Extracellular / circulating',
    'high-density lipoprotein particle': 'Extracellular / circulating',
    'triglyceride-rich plasma lipoprotein particle': 'Extracellular / circulating',
    'blood microparticle': 'Extracellular / circulating',
    'lipoprotein particle': 'Extracellular / circulating',
    'chylomicron': 'Extracellular / circulating',
    'collagen trimer': 'Extracellular / circulating',
    'focal adhesion': 'Extracellular / circulating',
    'cell-substrate junction': 'Extracellular / circulating',
    'integrin complex': 'Extracellular / circulating',
    'protein complex involved in cell adhesion': 'Extracellular / circulating',
    'adherens junction': 'Extracellular / circulating',
    'perisynaptic extracellular matrix': 'Extracellular / circulating',
    'synapse-associated extracellular matrix': 'Extracellular / circulating',
    # Intracellular
    'receptor complex': 'Intracellular',
    'plasma membrane raft': 'Intracellular',
    'membrane raft': 'Intracellular',
    'membrane microdomain': 'Intracellular',
    'plasma membrane signaling receptor complex': 'Intracellular',
    'catenin complex': 'Intracellular',
    'A band': 'Intracellular',
    'semaphorin receptor complex': 'Intracellular - degradation',
    'endocytic vesicle membrane': 'Intracellular - degradation',
    'endocytic vesicle': 'Intracellular - degradation',
    'early endosome': 'Intracellular - degradation',
    'endosome lumen': 'Intracellular - degradation',
    'clathrin-coated vesicle membrane': 'Intracellular - degradation',
    'clathrin-coated endocytic vesicle membrane': 'Intracellular - degradation',
    'clathrin-coated endocytic vesicle': 'Intracellular - degradation',
    'endolysosome': 'Intracellular - degradation',
    'endocytic vesicle lumen': 'Intracellular - degradation',
    'lytic vacuole': 'Intracellular - degradation',
    'vacuolar lumen': 'Intracellular - degradation',
    'lysosomal lumen': 'Intracellular - degradation',
    'peroxisomal matrix': 'Intracellular - degradation',
    'microbody lumen': 'Intracellular - degradation',
    'Proteasome regulatory particle, base subcomplex': 'Intracellular - degradation'
}

high_effect_cc_to_plot_sig['Theme'] = high_effect_cc_to_plot_sig['annotation_name'].map(theme_assignment)
high_effect_cc_to_plot_sig['Theme'] = high_effect_cc_to_plot_sig.Theme.fillna('Misc.')
high_effect_cc_to_plot_sig[
    'gene_ratio'] = high_effect_cc_to_plot_sig.n_test_annotated / high_effect_cc_to_plot_sig.n_total_annotated
high_effect_cc_to_plot_sig = high_effect_cc_to_plot_sig.sort_values(by='gene_ratio', ascending=False)
high_effect_cc_to_plot_sig['annotation_name'] = pd.Categorical(high_effect_cc_to_plot_sig['annotation_name'],
                                                               high_effect_cc_to_plot_sig['annotation_name'].to_list())

fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=go_c_fdr_sig_df, x='gene_ratio', y='annotation_name', hue='Theme', palette='tab20',
            hue_order=['Extracellular /\ncirculating', 'Secretory', 'Intracellular', 'Misc.'])
ax.set_ylabel('Gene Ontology (Cellular Component)')
ax.set_xlabel('Gene ratio')
ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
plt.tight_layout()
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/high_effect_olink_enriched_go_cc_terms.svg')
fig.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/high_effect_olink_enriched_go_cc_terms.png')
plt.close()
