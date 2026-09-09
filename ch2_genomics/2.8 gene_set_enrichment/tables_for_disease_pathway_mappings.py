""" Generate tables mapping diseases to Reactome terms corresponding to MAGMA FDR-sig enrichments. """
import pandas as pd

# background data
disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(str) + ')'

# results data
adjustment = 'minimal'

results = pd.read_csv(f"~/data/internal/genomics/magma/{adjustment}/reactome_enrichment_results_assoc.csv")
results = results[results.fdr_sig == 1][['disease', 'reactome_pathway_name']].merge(disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(columns='disease').rename(columns={'icd10_three_letter': 'disease'})

results_collated_list = []
for _, df in results.groupby('disease'):
    results_collated_list.append(df.assign(reactomes='; '.join(df.reactome_pathway_name.to_list())).drop(columns='reactome_pathway_name').drop_duplicates())
results_collated = pd.concat(results_collated_list)
results_collated[['disease', 'reactomes']].to_csv(f"~/data/internal/genomics/magma/{adjustment}/disease_reactome_mappings.csv", index=False)

