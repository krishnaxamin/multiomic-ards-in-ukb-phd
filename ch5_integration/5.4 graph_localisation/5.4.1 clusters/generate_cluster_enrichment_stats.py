""" Generate counts of clusters enriched in GWAS-only, proteomics-only, both, or neither. """
import pandas as pd
import os


def count_stats_per_summary(summary_df):

    summary_sig = summary_df[summary_df['sig_level'] == 'fdr'].copy()
    neither_count = summary_sig.neither_count.sum()
    gwas_count = summary_sig.gwas_count.sum()
    proteomics_count = summary_sig.proteomics_count.sum()
    both_count = summary_sig.both_count.sum()

    return {'neither': neither_count, 'gwas': gwas_count, 'proteomics': proteomics_count, 'both': both_count}


summaries_files = [x.name for x in
                   os.scandir(f"~/data/internal/knowledge_graph/gene_protein_distances/clustering") if
                   x.is_file() and 'summaries' in x.name]

neither_counter = 0
gwas_counter = 0
proteomics_counter = 0
both_counter = 0
for summary_file_name in summaries_files:

    summary_df = pd.read_csv(
            f"~/data/internal/knowledge_graph/gene_protein_distances/clustering/{summary_file_name}")
    counts = count_stats_per_summary(summary_df)

    neither_counter += counts['neither']
    gwas_counter += counts['gwas']
    proteomics_counter += counts['proteomics']
    both_counter += counts['both']

total_count = neither_counter + gwas_counter + proteomics_counter + both_counter

print(f"Total clusters: {total_count}\n"
      f"GWAS clusters: {gwas_counter} ({gwas_counter * 100 / total_count}%)\n"
      f"proteomics clusters: {proteomics_counter} ({proteomics_counter * 100 / total_count}%)\n"
      f"both clusters: {both_counter} ({both_counter * 100 / total_count}%)")

