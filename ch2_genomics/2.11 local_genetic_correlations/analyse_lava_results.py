import pandas as pd
import re
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import Dict, Union, List

from utils.significance_labelling import significance_labelling

mpl.use('TkAgg')

lava = pd.read_csv('~/data/internal/genomics/local_genetic_correlations/lava_ld_block_correlation.csv')

# remove MHC region LD blocks
ld_blocks_with_mhc = ['block' + str(x) for x in range(652, 659)]
lava = lava[~lava['locus'].isin(ld_blocks_with_mhc)].copy()

lava = significance_labelling(lava, fdr_permissive_group_by='locus')

lava_sig = lava[lava['fdr_sig'] == 1].copy()
# data descriptions
print(f"Number of LD blocks: {lava_sig.locus.nunique()}")
print(f"Number of diseases covered: {len(set(lava_sig.disease1) | set(lava_sig.disease2))}")

""" Plot scatterplot-heatmap """


# hue = rho, shape = significance (FDR = circle; FDR permissive = square; non-sig = triangle), size = n(shared loci)
# the more stringent FDR level is taken as preferred, i.e. if a disease pair is both FDR and FDR permissive, the former is taken
# (+) corr = upper half; (-) corr = lower half


def conditional_aggregate(group):
    """
    Averages over corr if all data points are the same sig_status. If there a mix of significance, the 'most sig' data
    points are averaged and the resulting average is deemed that highest sig_status.
    :param group:
    :return:
    """

    # if there is a mix, this picks the most significant status; otherwise, it gets the status all points have
    fdr_sig_status = group.fdr_sig.sort_values().iloc[0]
    fdr_sig_status_subset = group[group['fdr_sig'] == fdr_sig_status].copy()
    corr_avg = fdr_sig_status_subset['rho'].mean()
    n_values = len(fdr_sig_status_subset)

    return pd.Series({
        'mean_rho': corr_avg,
        'fdr_sig': fdr_sig_status,
        'n_blocks': n_values
    })


disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

plot_path = '~/ch2_genomics/2.11 local_genetic_correlations/plots/lava_correlation_lifestyles'

data_to_plot = lava[['disease1', 'disease2', 'rho', 'locus', 'fdr_sig']].copy()

data_to_plot['fdr_sig'] = data_to_plot['fdr_sig'].replace(1, 'Yes').replace(0, 'No')
data_to_plot = (data_to_plot
                .merge(disease_info[['disease_field', 'code_chapter']]
                       .rename(columns={'disease_field': 'disease1'}))
                .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                .merge(disease_info[['disease_field', 'code_chapter']]
                       .rename(columns={'disease_field': 'disease2'}))
                .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'})).drop_duplicates()
data_to_plot['fdr_sig'] = pd.Categorical(data_to_plot['fdr_sig'], categories=['Yes', 'No'])

# aggregate correlation and loci
data_to_plot_pos = data_to_plot[data_to_plot['rho'] > 0].copy()
data_to_plot_neg = data_to_plot[data_to_plot['rho'] < 0].rename(
    columns={'disease1': 'disease2', 'disease2': 'disease1'}).copy()
data_to_plot_pos_corr_agg = data_to_plot_pos.groupby(['disease1', 'disease2']).apply(
    conditional_aggregate).reset_index()
data_to_plot_neg_corr_agg = data_to_plot_neg.groupby(['disease1', 'disease2']).apply(
    conditional_aggregate).reset_index()

data_to_plot = pd.concat([data_to_plot_pos_corr_agg, data_to_plot_neg_corr_agg])
data_to_plot['disease1'] = pd.Categorical(data_to_plot['disease1'],
                                          categories=sorted(data_to_plot['disease1'].unique()))
data_to_plot['disease2'] = pd.Categorical(data_to_plot['disease2'],
                                          categories=sorted(data_to_plot['disease2'].unique()))

# map data to symmetric categorical index space
all_labels = sorted(set(data_to_plot['disease1']) | set(data_to_plot['disease2']))
# make mapping have 1-tick gaps between chapters, so that there are whitespace gaps/boundaries between chapters
label_to_num = {}
current_pos = 0
for i in range(len(all_labels)):
    if i == 0:
        label_to_num[all_labels[i]] = current_pos
        current_pos += 1
        continue
    current_chapter = int(re.search(r"\((\d+)\)", all_labels[i]).group(1))
    past_chapter = int(re.search(r"\((\d+)\)", all_labels[i - 1]).group(1))
    if past_chapter < current_chapter:
        current_pos += 1  # insert a 1-tick wide gap into the axis at the boundary between chapters
    label_to_num[all_labels[i]] = current_pos
    current_pos += 1
# label_to_num = {lab: i for i, lab in enumerate(all_labels)}
data_to_plot['disease1_as_num'] = data_to_plot['disease1'].map(label_to_num)
data_to_plot['disease2_as_num'] = data_to_plot['disease2'].map(label_to_num)

# plot
fig, ax = plt.subplots(figsize=(12, 10))
sns.scatterplot(data=data_to_plot.rename(
    columns={'fdr_sig': 'FDR < 0.05', 'mean_rho': 'Mean local correlation', 'n_blocks': 'Number LD blocks shared'}
),
    x='disease1_as_num', y='disease2_as_num', hue='Mean local correlation', style='FDR < 0.05',
    style_order=['Yes', 'No'], zorder=3, palette=plt.cm.coolwarm,
    hue_norm=(-1, 1), edgecolor=None, markers={'Yes': 'o', 'No': '^'},
    size='Number LD blocks shared', sizes=(20, 120), ax=ax)
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.ylabel('Disease ICD-10 code (ICD-10 chapter)', size=11)
plt.grid(True, zorder=0)

if ax.yaxis_inverted():
    ax.invert_yaxis()

# set ticks labels to the string versions
tick_positions = [label_to_num[label] for label in all_labels]
ax.set_xticks(tick_positions)
ax.set_xticklabels(all_labels, rotation=90, size=8)
ax.set_yticks(tick_positions)
ax.set_yticklabels(all_labels, size=8)

# add diagonal line going top-left to bottom-right, ensuring whitespaces are dealt with properly
coords = np.array(tick_positions)
ax.plot(coords, coords,
        color='black', linestyle='--', linewidth=0.5)

# bold legend section titles
legend = ax.get_legend()
for text, handle in zip(legend.texts, legend.legend_handles):
    if handle._label in ['Mean local correlation', 'Number LD blocks shared', 'FDR < 0.05']:  # section header
        text.set_fontweight('bold')

plt.tight_layout()
plt.savefig(f"{plot_path}.png")
plt.savefig(f"{plot_path}.svg")
plt.close()

""" Identify genes in the FDR-sig blocks"""


# try and do the FDR-sig-permissive genes, but if too many, fall back on FDR-sig only
# ignore above - doing only FDR-sig, for statistical rigour


def read_in_gene_set() -> pd.DataFrame:
    """
    Read in Ensembl's hg19 gene set.
    :return:
    """
    ensembl_hg19 = \
        pd.read_csv('https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz',
                    comment='#', sep='\t',
                    header=None,
                    names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                    usecols=['chr', 'feature', 'start', 'end', 'attribute'])[
            ['chr', 'feature', 'start', 'end', 'attribute']]
    ensembl_hg19 = ensembl_hg19[(~ensembl_hg19['chr'].isin(['Y', 'MT'])) & (ensembl_hg19['feature'] == 'gene')].drop(
        columns='feature').copy()

    # process attribute column to generate multiple columns for each attribute listed
    regex_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    attributes = ensembl_hg19['attribute'].map(
        lambda s: dict(regex_pattern.findall(str(s).strip().strip('"').replace('""', '"')))
    ).apply(pd.Series)
    ensembl_hg19 = pd.concat([ensembl_hg19, attributes], axis=1).drop(columns='attribute')

    return ensembl_hg19


def read_in_ld_blocks():
    """
    Read in Berisa and Pickrell 2016's LD blocks.
    :return:
    """

    ld_blocks = pd.read_csv('ukbiobank/gwas/berisa_pickrell_ld_blocks.bed', sep='\s+')
    ld_blocks['chr'] = ld_blocks['chr'].str.replace('chr', '')
    ld_blocks['locus'] = ['block' + str(x) for x in list(ld_blocks.index)]

    return ld_blocks


def map_ld_block_to_genes(block_id: str, ld_block_coord_rows: List[Dict[str, Union[int, str]]],
                          gene_set: pd.DataFrame) -> pd.Series:
    """
    Given an LD block from Berisa and Pickrell 2016, identify the Ensembl hg19 gene within it.
    Mapped genes must be fully within the LD block.
    :param block_id:
    :param ld_block_coord_rows:
    :param gene_set:
    :return:
    """

    block_coords = ld_block_coord_rows[int(block_id.replace('block', ''))]

    chromo = block_coords['chr']
    block_start = block_coords['start']
    block_end = block_coords['stop']

    mapped_genes = gene_set[(gene_set['chr'].astype(str) == chromo) & (gene_set['start'] >= block_start) & (
            gene_set['end'] <= block_end)].copy()
    mapped_genes_dict = {'ensg_id': '|'.join(mapped_genes.gene_id.to_list()),
                         'gene_symbol': '|'.join(mapped_genes.gene_name.to_list())}

    return pd.Series(mapped_genes_dict)


def map_lava_results_to_genes(
        lava_results_df: pd.DataFrame, ld_block_set: pd.DataFrame, gene_set: pd.DataFrame
) -> pd.DataFrame:
    """
    Given a set of LAVA results which contain block IDs in the 'locus' column, get Ensembl hg19 genes within those blocks.
    :param ld_block_set:
    :param lava_results_df:
    :param gene_set:
    :return:
    """
    ld_block_set_rows = ld_block_set.to_dict(orient='records')

    ld_blocks_to_map = lava_results_df[['locus']].drop_duplicates()
    ld_blocks_to_map[['mapped_ensg_ids', 'mapped_gene_symbols']] = ld_blocks_to_map['locus'].apply(
        map_ld_block_to_genes, ld_block_coord_rows=ld_block_set_rows, gene_set=gene_set)
    ld_blocks_to_map['mapped_ensg_ids'] = ld_blocks_to_map['mapped_ensg_ids'].replace('', 'Unmapped')
    ld_blocks_to_map['mapped_gene_symbols'] = ld_blocks_to_map['mapped_gene_symbols'].replace('', 'Unmapped')

    output = lava_results_df.merge(ld_blocks_to_map)

    return output


lava_sig = lava[lava['fdr_sig'] == 1].copy()
lava_sig_slim = lava_sig[['disease1', 'disease2', 'locus']].drop_duplicates().copy()

ensembl_hg19 = read_in_gene_set()
ld_blocks = read_in_ld_blocks()

lava_sig_mapped = map_lava_results_to_genes(lava_results_df=lava_sig_slim,
                                            ld_block_set=ld_blocks,
                                            gene_set=ensembl_hg19)

next_lava = lava_sig_mapped.copy()
next_lava['gene_id'] = next_lava['mapped_ensg_ids'].str.split('|')
next_lava = next_lava.explode('gene_id').drop(columns=['mapped_ensg_ids', 'mapped_gene_symbols'])
next_lava = next_lava.merge(
    ensembl_hg19[ensembl_hg19['gene_biotype'].isin(['protein_coding', 'lincRNA', 'antisense', 'miRNA'])]
    [['gene_id', 'chr', 'start', 'end', 'gene_biotype']])
next_lava.to_csv('~/ch2_genomics/2.11 local_genetic_correlations/int_data/loci_follow_up_genes_biotyped.csv', index=False)

genes_to_lava_loci = (next_lava[['gene_id', 'chr', 'start', 'end']]
                      .set_axis(['LOC', 'CHR', 'START', 'STOP'], axis=1)
                      .drop_duplicates())
genes_to_lava_loci.to_csv('~/ch2_genomics/2.11 local_genetic_correlations/int_data/loci_follow_up_genes.txt', sep='\t',
                          index=False)

gene_disease_pair = next_lava[['disease1', 'disease2', 'gene_id']].copy()
gene_disease_pair['disease_pair'] = gene_disease_pair['disease1'] + '-' + gene_disease_pair['disease2']
gene_disease_pair = gene_disease_pair.drop_duplicates()
(gene_disease_pair
 .to_csv('~/ch2_genomics/2.11 local_genetic_correlations/int_data/follow_up_genes_with_diseases.csv', index=False))

""" Import gene-level results """
lava_genes = pd.read_csv('~/data/internal/genomics/local_genetic_correlations/lava_gene_correlation.csv')
lava_genes = lava_genes[lava_genes['pval'].notna()].copy()

gene_disease_pair = pd.read_csv('~/ch2_genomics/2.11 local_genetic_correlations/int_data/follow_up_genes_with_diseases.csv')
next_lava = pd.read_csv('~/ch2_genomics/2.11 local_genetic_correlations/int_data/loci_follow_up_genes_biotyped.csv')
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
ensembl_hg19 = read_in_gene_set()

# remove disease pairs that are not relevant -
#  only want disease pair-gene combos where the gene was in an LD block for which there was FDR-sig correlation between the diseases
lava_genes = lava_genes.merge(gene_disease_pair.rename(columns={'gene_id': 'locus'})).sort_values(by='disease_pair')
lava_genes_with_blocks = lava_genes.merge(next_lava.rename(columns={'locus': 'block', 'gene_id': 'locus'}))

lava_genes_with_blocks = significance_labelling(lava_genes_with_blocks, fdr_permissive_group_by='locus')
lava_genes_sig = lava_genes_with_blocks[lava_genes_with_blocks['fdr_sig_permissive'] == 1][
    ['disease1', 'disease2', 'rho', 'locus']].copy()
lava_genes_sig = lava_genes_sig.merge(
    ensembl_hg19[['chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_biotype']].rename(columns={'gene_id': 'locus'}))
lava_genes_sig['gene_name'] = lava_genes_sig['gene_name'].fillna('Unnamed')
lava_genes_sig = lava_genes_sig.merge(next_lava.rename(columns={'locus': 'block', 'gene_id': 'locus'}))
lava_genes_sig = (lava_genes_sig
                  .merge(disease_info[['disease_field', 'code_chapter']]
                         .rename(columns={'disease_field': 'disease1'}))
                  .drop(columns='disease1').rename(columns={'code_chapter': 'disease1'})
                  .merge(disease_info[['disease_field', 'code_chapter']]
                         .rename(columns={'disease_field': 'disease2'}))
                  .drop(columns='disease2').rename(columns={'code_chapter': 'disease2'}))

# check heritability data - to look at why there are so few gene-wise correlations calculated
lava_heritability_genes = pd.read_csv(
    'ukbiobank/gwas/saige/genetic_similarities/lava/lava_gene_heritability_lifestyles.csv')
heritability_filtered = gene_disease_pair.drop(columns='disease_pair').melt(id_vars='gene_id',
                                                                            value_name='disease').drop(
    columns='variable').merge(lava_heritability_genes.rename(columns={'locus': 'gene_id'}))

# data descriptions
print(f"Number of genes: {lava_genes_sig.locus.nunique()}")
print(
    f"Breakdown of gene biotype: \n{lava_genes_sig[['locus', 'gene_biotype']].drop_duplicates().gene_biotype.value_counts()}")
print(f"Number of disease pairs: {(lava_genes_sig['disease1'] + lava_genes_sig['disease2']).nunique()}")

n_pairs_with_genes_sharing = lava_genes_sig.value_counts(['disease1', 'disease2'])
n_disease_pairs_per_gene = lava_genes_sig.value_counts(
    'gene_name')  # note there are 175 genes, but 174 of them have names

""" LAVA genes plotting """


def lava_gene_lollipop_plotting(lava_genes_sig_one_block, plot_path):
    """
    Plots spatially resolved LAVA local correlations (lollipop plot on top of gene tracks).
    Example outputs in gwas/saige/genetic_similarities/plots/lava/.
    :param lava_genes_sig_one_block:
    :param plot_path:
    :return:
    """

    lava_genes_sig_one_block['mid'] = (lava_genes_sig_one_block.start + lava_genes_sig_one_block.end) / 2
    lava_genes_sig_one_block['disease_pair'] = lava_genes_sig_one_block.disease1 + '-' + lava_genes_sig_one_block.disease2
    lava_genes_sig_one_block = lava_genes_sig_one_block.sort_values(by='disease_pair', ascending=False)

    pairs = sorted(list(lava_genes_sig_one_block.disease_pair.unique()), reverse=True)
    pair_index = {p: i for i, p in enumerate(pairs)}
    lava_genes_sig_one_block['y'] = lava_genes_sig_one_block.disease_pair.map(pair_index)

    genes = lava_genes_sig_one_block[['gene_name', 'start', 'end']].drop_duplicates().copy()
    genes = genes.sort_values('start')

    gene_cmap = plt.get_cmap('tab10')  # good categorical palette
    gene_colours = {
        gene: gene_cmap(i % gene_cmap.N)
        for i, gene in enumerate(genes['gene_name'])
    }

    fig, (ax1, ax2) = plt.subplots(
        2, 1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={'height_ratios': [3, 1.5]}
    )

    for _, r in lava_genes_sig_one_block.iterrows():
        gene_colour = gene_colours[r.gene_name]

        ax1.plot([r.mid, r.mid], [r.y, r.y + 0.4], lw=1.2, color=gene_colour)

        sc = ax1.scatter(
            r['mid'],
            r['y'] + 0.4,
            c=r['rho'],
            cmap='coolwarm',
            norm=mpl.colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
            s=60
        )

    ax1.set_yticks(range(len(pairs)))
    ax1.set_yticklabels(pairs)
    ax1.set_ylabel('Disease pairs: Disease ICD-10 code (ICD-10 chapter)')
    ax1.grid(
        axis='y',
        linestyle=':',
        linewidth=0.6,
        alpha=0.6
    )

    tracks = []
    track_ids = []
    for _, g in genes.iterrows():

        placed = False
        for i, track_end in enumerate(tracks):

            if g.start > track_end + 100000:
                tracks[i] = g.end
                track_ids.append(i)
                placed = True
                break
        if not placed:
            tracks.append(g.end)
            track_ids.append(len(tracks) - 1)

    genes['track'] = [x * 2 for x in track_ids]

    for _, g in genes.iterrows():
        y = g.track
        gene_colour = gene_colours[g.gene_name]

        ax2.plot(
            [g.start, g.end],
            [y, y],
            linewidth=8,
            color=gene_colour
        )

        ax2.text(
            (g.start + g.end) / 2,
            y + 1,
            g.gene_name,
            ha='center',
            va='center',
            fontsize=8
        )

    ax2.set_yticks([])
    ax2.set_xlabel('Genomic coordinate')
    ax2.set_ylabel('Genes')
    ax2.set_ylim([0, genes.track.max() + 2])

    plt.tight_layout()

    cbar = fig.colorbar(
        sc,
        ax=[ax1, ax2],
        location='right'
    )
    cbar.set_label('Genetic correlation')

    fig.savefig(f"{plot_path}.svg")
    fig.savefig(f"{plot_path}.png")

    plt.close()


for block in lava_genes_sig.block.unique():
    lava_gene_lollipop_plotting(lava_genes_sig_one_block=lava_genes_sig[lava_genes_sig['block'] == block].copy(),
                                plot_path=f"~/ch2_genomics/2.11 local_genetic_correlations/plots/lava/{block}")
