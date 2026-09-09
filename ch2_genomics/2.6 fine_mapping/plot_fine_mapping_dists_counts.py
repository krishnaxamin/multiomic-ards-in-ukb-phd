""" Generate plots and datasets (used to make tables) for fine mapping results. """
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import collections

mpl.use('TkAgg')

adjustment = 'lifestyles'

disease_info = pd.read_csv('ukbiobank/ard_identification/by_age_of_onset_stats/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(str) + ')'

processed_fine_mapping = pd.read_csv(
    f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/peak_pip0.5_fine_mapped_variants_blocksmerged_recredibled.csv")
processed_fine_mapping['mapped_gene_symbols'] = processed_fine_mapping['mapped_gene_symbols'].str.split('|')
processed_fine_mapping = processed_fine_mapping.explode('mapped_gene_symbols')

susie = processed_fine_mapping[processed_fine_mapping['algo'] == 'SuSiE'].copy()
abf = processed_fine_mapping[processed_fine_mapping['algo'] == 'ABF'].copy()

""" Number of regions/genes/variants per disease """
# ABF on LHS, SuSiE on RHS, medium Confidence top, high Confidence bottom
# for SuSiE, include all mapped vars, not just the top in the credible set

# ABF, medium
variants_abf_med = pd.DataFrame(abf[['disease', 'variant']].drop_duplicates().value_counts('disease')).reset_index()
genes_abf_med = pd.DataFrame(abf[abf.mapped_gene_symbols != 'Unmapped'][['disease', 'mapped_gene_symbols']].drop_duplicates().value_counts('disease')).reset_index()
regions_abf_med = pd.DataFrame(abf[['disease', 'merged_ld_region']].drop_duplicates().value_counts('disease')).reset_index()

# ABF, high
variants_abf_high = pd.DataFrame(abf[abf.pip > 0.9][['disease', 'variant']].drop_duplicates().value_counts('disease')).reset_index()
genes_abf_high = pd.DataFrame(abf[(abf.pip > 0.9) & (abf.mapped_gene_symbols != 'Unmapped')][['disease', 'mapped_gene_symbols']].drop_duplicates().value_counts('disease')).reset_index()
regions_abf_high = pd.DataFrame(abf[abf.pip > 0.9][['disease', 'merged_ld_region']].drop_duplicates().value_counts('disease')).reset_index()

# SuSiE, medium
variants_susie_med = pd.DataFrame(susie[['disease', 'variant']].drop_duplicates().value_counts('disease')).reset_index()
genes_susie_med = pd.DataFrame(susie[susie.mapped_gene_symbols != 'Unmapped'][['disease', 'mapped_gene_symbols']].drop_duplicates().value_counts('disease')).reset_index()
regions_susie_med = pd.DataFrame(susie[['disease', 'merged_ld_region']].drop_duplicates().value_counts('disease')).reset_index()

# SuSiE, high
variants_susie_high = pd.DataFrame(susie[susie.pip > 0.9][['disease', 'variant']].drop_duplicates().value_counts('disease')).reset_index()
genes_susie_high = pd.DataFrame(susie[(susie.pip > 0.9) & (susie.mapped_gene_symbols != 'Unmapped')][['disease', 'mapped_gene_symbols']].drop_duplicates().value_counts('disease')).reset_index()
regions_susie_high = pd.DataFrame(susie[susie.pip > 0.9][['disease', 'merged_ld_region']].drop_duplicates().value_counts('disease')).reset_index()

count_df = pd.concat([
    pd.concat([variants_abf_med.assign(Entity='Variant'),
               genes_abf_med.assign(Entity='Gene'),
               regions_abf_med.assign(Entity='Region')]).assign(Confidence='Medium', Algorithm='ABF'),
    pd.concat([variants_abf_high.assign(Entity='Variant'),
               genes_abf_high.assign(Entity='Gene'),
               regions_abf_high.assign(Entity='Region')]).assign(Confidence='High', Algorithm='ABF'),
    pd.concat([variants_susie_med.assign(Entity='Variant'),
               genes_susie_med.assign(Entity='Gene'),
               regions_susie_med.assign(Entity='Region')]).assign(Confidence='Medium', Algorithm='SuSiE'),
    pd.concat([variants_susie_high.assign(Entity='Variant'),
               genes_susie_high.assign(Entity='Gene'),
               regions_susie_high.assign(Entity='Region')]).assign(Confidence='High', Algorithm='SuSiE')
])

count_df = count_df.merge(disease_info[['icd10_three_letter', 'code_chapter']].set_axis(['disease', 'disease1'], axis=1)).drop(columns='disease').rename(columns={'disease1': 'disease', 'Confidence': 'Conf.'})
count_df['disease'] = pd.Categorical(count_df.disease, sorted(list(count_df.disease.unique())))

g = sns.catplot(data=count_df, x='count', y='disease', col='Algorithm', row='Conf.', kind='bar', sharey=True,
                sharex=True, hue='Entity', hue_order=['Variant', 'Gene', 'Region'], height=6, aspect=0.67)
g.set_axis_labels('Number of fine-mapped variants', 'Disease ICD-10 code (ICD-10 chapter)')
# g.set_titles("{col_name}", "{row_name}", size=15)
for ax in g.axes.flat:
    yticks = ax.get_yticks()
    # alternate shading
    for i in range(len(yticks)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, alpha=0.1, color='lightgrey')
    # add horizontal lines between categories
    for i in range(len(yticks) - 1):
        ax.axhline(i + 0.5, color='grey', linewidth=0.5, alpha=0.5)
    ax.set_ylim(-0.6, len(yticks) - 0.4)  # remove additional padding, while retaining a thin one
    ax.invert_yaxis()
    # add vertical lines on tick marks
    ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
    ax.set_axisbelow(True)
plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/n_vars_genes_regions_per_disease.svg")
plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/n_vars_genes_regions_per_disease.png")
plt.close()

""" Generate disease-gene table """
genes_only = processed_fine_mapping[['disease', 'pip', 'algo', 'mapped_gene_symbols']].copy()
genes_only = genes_only[genes_only.mapped_gene_symbols != 'Unmapped'].copy()
genes_only = genes_only.sort_values(by=['disease', 'algo', 'pip'], ascending=[True, True, False]).drop_duplicates(['disease', 'algo', 'mapped_gene_symbols'])
genes_only.loc[genes_only.pip > 0.9, 'pip'] = 'High'
genes_only.loc[genes_only.pip != 'High', 'pip'] = 'Medium'

genes_only_collated_list = []
for _, df in genes_only.groupby(['disease', 'algo', 'pip']):
    genes_only_collated_list.append(df.assign(genes=', '.join(df.mapped_gene_symbols.to_list())).drop(columns='mapped_gene_symbols').drop_duplicates())
genes_only_collated = pd.concat(genes_only_collated_list)
genes_only_collated = genes_only_collated[['disease', 'pip', 'genes', 'algo']].rename(columns={'pip': 'confidence'})
genes_only_collated.to_csv(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/disease_genes_mappings.csv", index=False)

""" Genes mapping to multiple diseases """
genes_only_code_chapter = genes_only.merge(disease_info[['icd10_three_letter', 'code_chapter']].rename(columns={'icd10_three_letter': 'disease'}))
genes_only_code_chapter['pip'] = genes_only_code_chapter['pip'].map({'High': 1, 'Medium': -1})
# heatmaps: genes on y, ARDS on x, ABF on LHS, SuSiE on RHS
abf_gene_freq = genes_only_code_chapter[genes_only_code_chapter.algo == 'ABF'].value_counts('mapped_gene_symbols')
abf_gene_freq_heatmap = genes_only_code_chapter[(genes_only_code_chapter.algo == 'ABF') &
                                                (genes_only_code_chapter.mapped_gene_symbols.isin(
                                                    list(abf_gene_freq[abf_gene_freq > 1].index)))].pivot(index='mapped_gene_symbols', columns='code_chapter', values='pip').fillna(0)

susie_gene_freq = genes_only_code_chapter[genes_only_code_chapter.algo == 'SuSiE'].value_counts('mapped_gene_symbols')
susie_gene_freq_heatmap = genes_only_code_chapter[(genes_only_code_chapter.algo == 'SuSiE') &
                                                (genes_only_code_chapter.mapped_gene_symbols.isin(
                                                    list(susie_gene_freq[susie_gene_freq > 1].index)))].pivot(index='mapped_gene_symbols', columns='code_chapter', values='pip').fillna(0)

joint_genes = sorted(list(set(abf_gene_freq_heatmap.index) | set(susie_gene_freq_heatmap.index)))
joint_diseases = sorted(list(set(abf_gene_freq_heatmap.columns) | set(susie_gene_freq_heatmap.columns)))

abf_gene_freq_heatmap = abf_gene_freq_heatmap.reindex(index=joint_genes, columns=joint_diseases, fill_value=0)
susie_gene_freq_heatmap = susie_gene_freq_heatmap.reindex(index=joint_genes, columns=joint_diseases, fill_value=0)

fig = plt.figure(figsize=(8, 6))
fig.set_constrained_layout(True)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)

cmap = mpl.colors.LinearSegmentedColormap.from_list('coolwarm_white', [plt.cm.coolwarm(0.0), 'white', plt.cm.coolwarm(1.0)])

ax = fig.add_subplot(gs[0])
sns.heatmap(
    data=abf_gene_freq_heatmap, linewidths=0.35, ax=ax, xticklabels=True, yticklabels=True,
    cbar=False, cmap=cmap, linecolor='black'
)
ax.set_title('ABF')
ax.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax.set_ylabel('Gene')
ax.set_xticklabels(ax.get_xticklabels(), size=9)
ax.set_yticklabels(ax.get_yticklabels(), size=9)

ax2 = fig.add_subplot(gs[1], sharey=ax)
sns.heatmap(
    data=susie_gene_freq_heatmap, linewidths=0.35, ax=ax2, xticklabels=True, yticklabels=True,
    cbar=False, cmap=cmap, linecolor='black'
)
ax2.set_title('SuSiE')
ax2.set_xlabel('Disease ICD-10 code (ICD-10 chapter)')
ax2.set_ylabel('')
ax2.set_xticklabels(ax.get_xticklabels(), size=9)
ax2.tick_params(axis='y', left=False, labelleft=False)  # hide y-axis labels on RHS plot

cmap = plt.get_cmap("coolwarm")
patches = [
    mpl.patches.Patch(color=cmap(1.0), label='High'),
    mpl.patches.Patch(color=cmap(0.0), label='Medium')
]
ax2.legend(handles=patches, bbox_to_anchor=(1.05, 0.5), loc='center left', title='Confidence')

plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/genes_with_multiple_diseases.png")
plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/genes_with_multiple_diseases.svg")
plt.close()

""" Compare SuSiE and ABF results """
# Sharing of blocks that have fine-mapped vars (e.g. one has 5 blocks, one has 3 blocks, 2 blocks are shared);
# sharing of genes (numbers and which ones are shared, which ones is susie specific);
# sharing of variants (numbers and which ones are shared, which ones are susie specific)
mapped_vars_with_merged_regions_recredibled = pd.read_csv(
    f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/peak_pip0.5_fine_mapped_variants_blocksmerged_recredibled.csv")
confidence = 'medium'

if confidence == 'medium':
    susie = mapped_vars_with_merged_regions_recredibled[
        (mapped_vars_with_merged_regions_recredibled['algo'] == 'SuSiE') & (
                    mapped_vars_with_merged_regions_recredibled['pip'] <= 0.9)].copy()
    abf = mapped_vars_with_merged_regions_recredibled[
        (mapped_vars_with_merged_regions_recredibled['algo'] == 'ABF') & (
                    mapped_vars_with_merged_regions_recredibled['pip'] <= 0.9)].copy()
elif confidence == 'high':
    susie = mapped_vars_with_merged_regions_recredibled[
        (mapped_vars_with_merged_regions_recredibled['algo'] == 'SuSiE') & (
                    mapped_vars_with_merged_regions_recredibled['pip'] > 0.9)].copy()
    abf = mapped_vars_with_merged_regions_recredibled[(mapped_vars_with_merged_regions_recredibled['algo'] == 'ABF') & (
                mapped_vars_with_merged_regions_recredibled['pip'] > 0.9)].copy()

susie_blocks = {disease: set(df['merged_ld_region']) for disease, df in susie.groupby('disease')}
abf_blocks = {disease: set(df['merged_ld_region']) for disease, df in abf.groupby('disease')}
shared_diseases = sorted(list(set(susie_blocks.keys()) & set(abf_blocks.keys())))
sharing_dict = collections.defaultdict(dict)

# identify which blocks/genes/variants are shared/unique + export
for disease in shared_diseases:
    shared_blocks = susie_blocks[disease] & abf_blocks[disease]
    n_shared_blocks = len(shared_blocks)
    only_susie_blocks = susie_blocks[disease] - abf_blocks[disease]
    n_only_susie_blocks = len(only_susie_blocks)
    only_abf_blocks = abf_blocks[disease] - susie_blocks[disease]
    n_only_abf_blocks = len(only_abf_blocks)

    sharing_dict[disease] = {'shared_blocks': '|'.join(shared_blocks),
                             'susie_only_blocks': '|'.join(only_susie_blocks),
                             'abf_only_blocks': '|'.join(only_abf_blocks),
                             'n_shared_blocks': n_shared_blocks,
                             'n_susie_only_blocks': n_only_susie_blocks,
                             'n_abf_only_blocks': n_only_abf_blocks}

    if n_only_abf_blocks + n_only_susie_blocks > 0:
        print(
            f"{disease}: n_shared_blocks={n_shared_blocks}; n_susie_only={n_only_susie_blocks}; n_abf_only={n_only_abf_blocks}")

susie_genes = {disease: set(df.mapped_gene_symbols.str.split('|').explode('mapped_gene_symbols')) - {'Unmapped'} for
               disease, df in susie.groupby('disease')}
abf_genes = {disease: set(df.mapped_gene_symbols.str.split('|').explode('mapped_gene_symbols')) - {'Unmapped'} for
             disease, df in abf.groupby('disease')}
print('Gene-level')

for disease in shared_diseases:
    shared_genes = susie_genes[disease] & abf_genes[disease]
    n_shared_genes = len(shared_genes)
    only_susie_genes = susie_genes[disease] - abf_genes[disease]
    n_only_susie_genes = len(only_susie_genes)
    only_abf_genes = abf_genes[disease] - susie_genes[disease]
    n_only_abf_genes = len(only_abf_genes)

    sharing_dict[disease].update({'shared_genes': '|'.join(shared_genes),
                                  'susie_only_genes': '|'.join(only_susie_genes),
                                  'abf_only_genes': '|'.join(only_abf_genes),
                                  'n_shared_genes': n_shared_genes,
                                  'n_susie_only_genes': n_only_susie_genes,
                                  'n_abf_only_genes': n_only_abf_genes})

    if n_only_abf_genes + n_only_susie_genes > 0:
        print(
            f"{disease}: n_shared_genes={n_shared_genes}; n_susie_only={n_only_susie_genes}; n_abf_only={n_only_abf_genes}")

susie_vars = {disease: set(df['variant']) for disease, df in susie.groupby('disease')}
abf_vars = {disease: set(df['variant']) for disease, df in abf.groupby('disease')}
print('Variant-level')

for disease in shared_diseases:
    shared_vars = susie_vars[disease] & abf_vars[disease]
    n_shared_vars = len(shared_vars)
    only_susie_vars = susie_vars[disease] - abf_vars[disease]
    n_only_susie_vars = len(only_susie_vars)
    only_abf_vars = abf_vars[disease] - susie_vars[disease]
    n_only_abf_vars = len(only_abf_vars)

    sharing_dict[disease].update({'shared_vars': '|'.join(shared_vars),
                                  'susie_only_vars': '|'.join(only_susie_vars),
                                  'abf_only_vars': '|'.join(only_abf_vars),
                                  'n_shared_vars': n_shared_vars,
                                  'n_susie_only_vars': n_only_susie_vars,
                                  'n_abf_only_vars': n_only_abf_vars})

    if n_only_abf_vars + n_only_susie_vars > 0:
        print(
            f"{disease}: n_shared_vars={n_shared_vars}; n_susie_only={n_only_susie_vars}; n_abf_only={n_only_abf_vars}")

sharing_dict_as_list = [{'disease': disease} | dictt for disease, dictt in sharing_dict.items()]
sharing_df = pd.DataFrame(sharing_dict_as_list)
sharing_df.to_csv(
    f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/susie_abf_concordance_counting_{confidence}_confidence.csv",
    index=False)

# plot sharing_df: just plot number of variants/genes/regions mapped by SuSiE only/ABF only/both
# have lines demarcating increments or text above each column
# can facet to plot both high and medium confidence mappings in one plot, as well as variants/genes/regions
# try bar chart, rather than col chart
disease_info = pd.read_csv(
    'ukbiobank/ard_identification/by_age_of_onset_stats/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'
sharing_df_to_plot = pd.concat([pd.read_csv(
    f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/susie_abf_concordance_counting_medium_confidence.csv").assign(
    conf_level='Medium'),
                                pd.read_csv(
                                    f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/susie_abf_concordance_counting_high_confidence.csv").assign(
                                    conf_level='High')])
sharing_df_to_plot = sharing_df_to_plot[
    ['disease', 'conf_level'] + [x for x in sharing_df_to_plot.columns if 'n_' in x]]
sharing_df_to_plot_long = pd.concat([
    sharing_df_to_plot[['disease', 'conf_level'] + [x for x in sharing_df_to_plot.columns if 'blocks' in x]].melt(
        id_vars=['disease', 'conf_level'], value_name='count', var_name='class').assign(entity='Region'),
    sharing_df_to_plot[['disease', 'conf_level'] + [x for x in sharing_df_to_plot.columns if 'genes' in x]].melt(
        id_vars=['disease', 'conf_level'], value_name='count', var_name='class').assign(entity='Gene'),
    sharing_df_to_plot[['disease', 'conf_level'] + [x for x in sharing_df_to_plot.columns if 'vars' in x]].melt(
        id_vars=['disease', 'conf_level'], value_name='count', var_name='class').assign(entity='Variant')
])
sharing_df_to_plot_long['class'] = ['_'.join(x.split('_')[:-1]) for x in list(sharing_df_to_plot_long['class'])]
sharing_df_to_plot_long['class'] = (sharing_df_to_plot_long['class']
                                    .replace(to_replace='n_shared', value='Both ABF and SuSiE')
                                    .replace(to_replace='n_susie_only', value='SuSiE only')
                                    .replace(to_replace='n_abf_only', value='ABF only'))
sharing_df_to_plot_long = sharing_df_to_plot_long.merge(
    disease_info[['icd10_three_letter', 'code_chapter']].rename(columns={'icd10_three_letter': 'disease'}))
sharing_df_to_plot_long['code_chapter'] = pd.Categorical(sharing_df_to_plot_long.code_chapter, sorted(list(set(sharing_df_to_plot_long.code_chapter))))

g = sns.catplot(
    data=sharing_df_to_plot_long.rename(columns={'conf_level': 'Conf.', 'entity': 'Entity', 'class': 'Algorithm'}),
    x='count', y='code_chapter', hue='Algorithm', col='Entity', row='Conf.', kind='bar', sharey=True,
    sharex=True, hue_order=['ABF only', 'SuSiE only', 'Both ABF and SuSiE'],
    height=6, aspect=0.67)
g.set_axis_labels('Number of fine-mapped variants', 'Disease ICD-10 code (ICD-10 chapter)')
# g.set_titles("{col_name}", "{row_name}", size=15)
for ax in g.axes.flat:
    yticks = ax.get_yticks()
    # alternate shading
    for i in range(len(yticks)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, alpha=0.1, color='lightgrey')
    # add horizontal lines between categories
    for i in range(len(yticks) - 1):
        ax.axhline(i + 0.5, color='grey', linewidth=0.5, alpha=0.5)
    ax.set_ylim(-0.6, len(yticks) - 0.4)  # remove additional padding, while retaining a thin one
    ax.invert_yaxis()
    # add vertical lines on tick marks
    ax.xaxis.grid(True, which='major', color='grey', alpha=0.2)
    ax.set_axisbelow(True)
plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/susie_abf_concordance.svg")
plt.savefig(f"ukbiobank/gwas/saige/fine_mapping/{adjustment}/plots/susie_abf_concordance.png")
plt.close()
