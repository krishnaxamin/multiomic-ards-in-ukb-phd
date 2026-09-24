""" Determine in which GTEx tissues genes have been identified as tissue-specific by the HPA. """
from pandas import read_csv, concat

expression = read_csv('~/data/external/human_protein_atlas/rna_tissue_detail_gtex.tsv', sep='\t')
tissue_specificity = read_csv('~/data/external/human_protein_atlas/tissue_category_rna_any_Tissue.tsv', sep='\t')

# 10926 genes display some sort of specificity
# Tissue enriched = >=4x expression vs next best
# Group enriched = >=4x expression in 2-5 tissues vs next best
# Tissue enhanced = >=4x expression vs average of all other tissues

""" Tissue enriched """
# Sort by class and mark descending
expression_sorted = expression.sort_values(['Gene', 'nTPM'], ascending=[True, False])

# For each gene, get top 2 tissues
top2 = expression_sorted.groupby('Gene').head(2)

# Keep only genes that actually have ≥ 2 tissues
valid_classes = top2['Gene'].value_counts()[lambda x: x >= 2].index
top2 = top2[top2['Gene'].isin(valid_classes)]

# Compute ratio top / second
ratios = top2.groupby('Gene')['nTPM'].apply(lambda x: x.iloc[0] / x.iloc[1])

# Find genes where top tissue >= 4 × second best
classes_extreme = ratios[ratios >= 4].index

# Extract genes in those tissues whose mark is the top by 4x
extremes = expression_sorted.groupby('Gene').head(1)
tissue_enriched = extremes[extremes['Gene'].isin(classes_extreme)]

tissue_enriched_annot = tissue_enriched.merge(tissue_specificity[tissue_specificity['RNA tissue specificity'] == 'Tissue enriched'].rename(columns={'Gene': 'Gene name'}))

""" Tissue enhanced """
# Compute gene averages
gene_avg = expression.groupby('Gene')['nTPM'].mean().rename('gene_avg')

# Merge averages back into expression
expression_with_avg = expression.merge(gene_avg, on='Gene')

# Sort by gene and nTPM descending
expression_sorted = expression_with_avg.sort_values(['Gene', 'nTPM'], ascending=[True, False])

# Take the top tissue per gene
top_tissues = expression_sorted.groupby('Gene').head(1)

# Compute ratio top / gene average
top_tissues['ratio'] = top_tissues['nTPM'] / top_tissues['gene_avg']

# Keep only those with ≥ 4× gene average
tissue_enhanced = top_tissues[top_tissues['ratio'] >= 4].drop(columns='gene_avg')

tissue_enhanced_annot = tissue_enhanced.merge(tissue_specificity[tissue_specificity['RNA tissue specificity'] == 'Tissue enhanced'].rename(columns={'Gene': 'Gene name'}))

""" Group enriched """


# For each Gene, identify 2-5 tissues where expression is >=4x higher than any other tissue
def find_extreme_tissues(sub):
    sub = sub.sort_values("nTPM", ascending=False)
    for k in range(2, 6):  # check top 2–5
        if len(sub) >= k + 1:  # need at least k+1 tissues to compare
            top = sub.head(k)
            rest = sub.iloc[k:]
            if (rest.empty) or (top["nTPM"].min() >= 4 * rest["nTPM"].max()):
                # First valid k → smallest valid set
                return top.assign(n_top=k)
    return None


# Sort tissues by nTPM within each Gene
expression_sorted = expression.sort_values(["Gene", "nTPM"], ascending=[True, False])

group_enriched = (
    expression_sorted.groupby("Gene", group_keys=False)
    .apply(find_extreme_tissues)
)

group_enriched_annot = group_enriched.merge(tissue_specificity[tissue_specificity['RNA tissue specificity'] == 'Group enriched'].rename(columns={'Gene': 'Gene name'}))

""" Combine all and export """
tissue_specificity_with_tissue_id = concat([tissue_enriched_annot,
                                            group_enriched_annot.drop(columns='n_top'),
                                            tissue_enhanced_annot.drop(columns='ratio')])
# 6826 genes can have their specificity mapped to a GTEx tissue
tissue_specificity_with_tissue_id.to_csv('~/data/external/human_protein_atlas/tissue_specificity_with_tissue_id.csv', index=False)
