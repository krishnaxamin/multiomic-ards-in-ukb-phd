""" SAIGE vs REGENIE comparisons """
import pandas as pd
import numpy as np
import scipy
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('TkAgg')

# read in disease info
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/genomics_qcv1_qu10_50_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Number of variants """

saige = pd.concat([pd.read_csv(
    '~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/pan-ukbb-eur_assoc_geno_saige.csv').assign(info_label='high_info'),
                   pd.read_csv(
                       '~/ch2_genomics/2.2 trial_gwas_runs/2.2.1 saige/pan-ukbb-eur_assoc_imputed_saige.csv')]).drop_duplicates(
    ['CHR', 'POS', 'REF', 'ALT', 'disease'])
regenie = pd.concat(
    [pd.read_csv('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/pan-ukbb-eur_assoc_geno_regenie.csv').assign(info_label='high_info'),
     pd.read_csv('~/ch2_genomics/2.2 trial_gwas_runs/2.2.2 regenie/pan-ukbb-eur_assoc_imputed_regenie.csv')]).drop_duplicates(
    ['CHROM', 'GENPOS', 'REF', 'ALT', 'disease'])

saige_assoc = saige[(saige['assoc_label'] == 'assoc') & (saige['info_label'] == 'high_info')].copy()
regenie_assoc = regenie[(regenie['LOG10P'] > -1 * np.log10(5e-8 / 69)) & (regenie['info_label'] == 'high_info')].copy()

# number of diseases with variants significantly associated
saige_assoc.disease.nunique()
regenie_assoc.disease.nunique()
len(set(saige_assoc.disease) & set(regenie_assoc.disease))

# number of variants significantly associated with each disease
vars_per_disease = pd.concat([saige_assoc.disease.value_counts(), regenie_assoc.disease.value_counts()], axis=1).set_axis(['saige', 'regenie'], axis=1).reset_index()
scipy.stats.wilcoxon(vars_per_disease.saige, vars_per_disease.regenie)
scipy.stats.wilcoxon(vars_per_disease.saige, vars_per_disease.regenie, alternative='less')
vars_per_disease['saige_to_regenie_gained'] = vars_per_disease['regenie'] - vars_per_disease['saige']
vars_per_disease['saige_to_regenie_gained_as_perc'] = vars_per_disease['saige_to_regenie_gained'] / vars_per_disease['saige']
vars_per_disease['saige_to_regenie_gained_as_perc'].describe()
(vars_per_disease['saige_to_regenie_gained'] > 0).sum()
(vars_per_disease['saige_to_regenie_gained'] < 0).sum()
(vars_per_disease['saige_to_regenie_gained'] == 0).sum()
vars_per_disease.to_csv('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/num_assoc_vars_per_disease.csv', index=False)

# number of unique variants
len(saige_assoc[['CHR', 'POS', 'Allele1', 'Allele2']].drop_duplicates())  # 29572
len(regenie_assoc[['CHROM', 'GENPOS', 'REF', 'ALT']].drop_duplicates())  # 30049
n_unique_vars = pd.concat([saige_assoc[['CHR', 'POS', 'Allele1', 'Allele2']].dri(),
                           regenie_assoc.disease[['CHROM', 'GENPOS', 'REF', 'ALT']].value_counts()], axis=1).set_axis(['saige', 'regenie'], axis=1).reset_index()

# plot
vars_per_disease_to_plot = vars_per_disease.merge(disease_info[['icd10_three_letter', 'code_chapter']].rename(columns={'icd10_three_letter': 'disease'})).drop(columns='disease').sort_values(by='saige_to_regenie_gained_as_perc', ascending=False)
vars_per_disease_to_plot['code_chapter'] = pd.Categorical(vars_per_disease_to_plot['code_chapter'], list(vars_per_disease_to_plot['code_chapter']))

plt.figure(figsize=(12, 8))
sns.barplot(x=vars_per_disease_to_plot['code_chapter'], y=vars_per_disease_to_plot['saige_to_regenie_gained_as_perc'] * 100)
plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
plt.xticks(rotation=90)
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('(nVars(REGENIE) - nVars(SAIGE)) / nVars(SAIGE), as percentage')
plt.tight_layout()
plt.savefig('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/plots/num_assoc_vars_per_disease.svg')
plt.savefig('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/plots/num_assoc_vars_per_disease.png')
plt.close()

""" GIFs """
gifs = pd.concat([pd.read_csv(f"~/ch2_genomics/2.3 genomic_inflation_factors/disease_gifs/{file}").assign(disease=file.split('_')[0]) for file in [x.name for x in os.scandir('~/ch2_genomics/2.3 genomic_inflation_factors/disease_gifs') if x.is_file()]])
scipy.stats.wilcoxon(gifs.saige_gif_pval, gifs.regenie_gif_pval)
scipy.stats.wilcoxon(gifs.saige_gif_pval, gifs.regenie_gif_pval, alternative='less')
(gifs.saige_gif_pval < gifs.regenie_gif_pval).sum()
gifs['saige_to_regenie_as_perc_of_saige_diff_from_1'] = (gifs.regenie_gif_pval - gifs.saige_gif_pval) / (gifs.saige_gif_pval - 1)
gifs['saige_to_regenie_as_perc_of_saige_diff_from_1'].describe()

# plot
gifs = gifs.rename(columns={'disease': 'disease_field'}).merge(disease_info[['disease_field', 'code_chapter']].rename(columns={'code_chapter': 'disease'})).sort_values(by='saige_to_regenie_as_perc_of_saige_diff_from_1', ascending=False)
gifs['disease'] = pd.Categorical(gifs['disease'], list(gifs['disease']))

plt.figure(figsize=(12, 8))
sns.barplot(x=gifs['disease'], y=gifs['saige_to_regenie_as_perc_of_saige_diff_from_1'] * 100)
plt.xticks(rotation=90)
plt.xlabel('Disease ICD-10 code (ICD-10 chapter)')
plt.ylabel('(GIF(REGENIE) - GIF(SAIGE)) / (GIF(SAIGE) - 1), as percentage')
plt.tight_layout()
plt.savefig('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/plots/saige_regenie_gifs.svg')
plt.savefig('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/plots/saige_regenie_gifs.png')
plt.close()

gifs.to_csv('~/ch2_genomics/2.2 trial_gwas_runs/saige_vs_regenie/saige_regenie_gifs.csv', index=False)

