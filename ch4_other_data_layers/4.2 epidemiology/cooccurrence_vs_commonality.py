""" Compare disease co-occurrence with various disease commonality results from different omics and cohorts. """
import pandas as pd
import statsmodels.api as sm

from utils.significance_labelling import significance_labelling

""" Regressions """


def linear(indep_var_df: pd.DataFrame, dep_var_df: pd.DataFrame, indep_var_col: str, dep_var_col: str):
    # make sure 1:1 between the dfs
    indep_var_df['disease_pair'] = indep_var_df.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                      axis=1)
    dep_var_df['disease_pair'] = dep_var_df.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                  axis=1)
    vars_df = indep_var_df[['disease_pair', indep_var_col]].merge(dep_var_df[['disease_pair', dep_var_col]])

    # fit model
    results = sm.OLS(vars_df[dep_var_col], vars_df[indep_var_col]).fit()

    return {'coeff': results.params[indep_var_col], 'pvalue': results.pvalues[indep_var_col], 'algo': 'linear'}


def logistic(indep_var_df: pd.DataFrame, dep_var_df: pd.DataFrame, indep_var_col: str, dep_var_col: str):
    # make sure all independent samples have a matched dependent outcome
    indep_var_df['disease_pair'] = indep_var_df.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                      axis=1)
    dep_var_df['disease_pair'] = dep_var_df.apply(lambda row: '-'.join(sorted([row['disease1'], row['disease2']])),
                                                  axis=1)
    vars_df = indep_var_df[['disease_pair', indep_var_col]].merge(dep_var_df[['disease_pair', dep_var_col]])

    # fit model
    results = sm.Logit(vars_df[dep_var_col], vars_df[indep_var_col]).fit()

    return {'coeff': results.params[indep_var_col], 'pvalue': results.pvalues[indep_var_col], 'algo': 'logistic'}


""" Invariant info """
disease_info = pd.read_csv(
    '~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
    str) + ')'

""" Genomics """
# alpha-MLE vs significant entity sharing & SemSim

# alpha-MLE
genomics = pd.read_csv('~/data/internal/genomics/genomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
genomics = significance_labelling(genomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                  pair_columns=['disease1', 'disease2'])
genomics_sig = genomics[genomics['fdr_sig'] == 1].copy()
genomics_sig = genomics_sig.loc[~genomics_sig.index.duplicated(keep='first'), :]

# variant sharing
variant_sharing = pd.read_csv('~/data/internal/genomics/shared_assoc_entities/variant_sharing_fisher_cossims.csv')
variant_sharing = (variant_sharing.merge(
    disease_info[['disease_field', 'icd10_three_letter']]
    .rename(columns={'icd10_three_letter': 'disease1'}))
                   .drop(columns='disease1').rename(columns={'disease_field': 'disease1'})
                   .merge(disease_info[['disease_field', 'icd10_three_letter']]
                          .rename(columns={'icd10_three_letter': 'disease2'}))
                   .drop(columns='disease2').rename(columns={'disease_field': 'disease2'}))

# gene sharing
gene_sharing = pd.read_csv('~/data/internal/genomics/shared_assoc_entities/gene_sharing_fisher.csv')

# HDL
hdl = pd.read_csv('~/data/internal/genomics/whole_genome_genetic_correlations/hdl_genetic_correlations.csv')
hdl = significance_labelling(hdl, pvalue_col='pval')

# SemSim
genomics_semsim = pd.read_csv('data/internal/genomics/pathway_based_similarities/disease_resnik_bma_reactome.csv')

# alpha vs sig-sharing-or-not
genomics_sig_variant_sharing_vs_alpha = logistic(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                                 dep_var_df=variant_sharing, dep_var_col='fdr_sig')
genomics_sig_gene_sharing_vs_alpha = logistic(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                              dep_var_df=gene_sharing, dep_var_col='fdr_sig')
# alpha vs number of shared entities
genomics_n_variant_sharing_vs_alpha = linear(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                             dep_var_df=variant_sharing[variant_sharing['fdr_sig'] == 1].copy(),
                                             dep_var_col='n_shared')
genomics_n_gene_sharing_vs_alpha = linear(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                          dep_var_df=gene_sharing[gene_sharing['fdr_sig'] == 1].copy(),
                                          dep_var_col='n_shared')
# alpha vs sig-HDL-or-not
genomics_sig_hdl_vs_alpha = logistic(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                     dep_var_df=hdl, dep_var_col='fdr_sig')
# alpha vs HDL
genomics_hdl_vs_alpha = linear(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                               dep_var_df=hdl[hdl['fdr_sig'] == 1].copy(),
                               dep_var_col='genetic_correlation')
# alpha vs SemSim
genomics_semsim_vs_alpha = linear(indep_var_df=genomics_sig, indep_var_col='alpha_hat',
                                  dep_var_df=genomics_semsim, dep_var_col='bma_sem_sim')
# store results
genomics_results = {'sig_variant_sharing': genomics_sig_variant_sharing_vs_alpha,
                    'sig_gene_sharing': genomics_sig_gene_sharing_vs_alpha,
                    'n_variant_sharing': genomics_n_variant_sharing_vs_alpha,
                    'n_gene_sharing': genomics_n_gene_sharing_vs_alpha,
                    'sig_hdl': genomics_sig_hdl_vs_alpha,
                    'hdl': genomics_hdl_vs_alpha,
                    'semsim': genomics_semsim_vs_alpha}

""" Prevalent/prior Proteomics """
# prevalent alpha-MLE
prior_proteomics = pd.read_csv(
    '~/data/internal/proteomics/proteomics_pan_ukbb_eur_prior_qu10_50_cooccurrences_alphamle.csv')
prior_proteomics = significance_labelling(prior_proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                          pair_columns=['disease1', 'disease2'])
prior_proteomics_sig = prior_proteomics[prior_proteomics['fdr_sig'] == 1].copy()
prior_proteomics_sig = prior_proteomics_sig.loc[~prior_proteomics_sig.index.duplicated(keep='first'), :]

# protein sharing
firth_protein_sharing = pd.read_csv(
    '~/data/internal/proteomics/shared_assoc_proteins/firth_protein_sharing_fisher_cossims.csv')
firth_protein_sharing = (firth_protein_sharing.merge(
    disease_info[['disease_field', 'code_chapter']]
    .rename(columns={'code_chapter': 'disease1'}))
                         .drop(columns='disease1').rename(columns={'disease_field': 'disease1'})
                         .merge(disease_info[['disease_field', 'code_chapter']]
                                .rename(columns={'code_chapter': 'disease2'}))
                         .drop(columns='disease2').rename(columns={'disease_field': 'disease2'}))

# SemSim
firth_semsim = pd.read_csv(
    '~/data/internal/proteomics/pathway_based_similarities/firth_disease_resnik_bma.csv')

# alpha vs sig-sharing-or-not
firth_sig_protein_sharing_vs_prevalent_alpha = logistic(indep_var_df=prior_proteomics_sig, indep_var_col='alpha_hat',
                                                        dep_var_df=firth_protein_sharing, dep_var_col='fdr_sig')
# alpha vs number of shared entities
firth_n_protein_sharing_vs_prevalent_alpha = linear(indep_var_df=prior_proteomics_sig, indep_var_col='alpha_hat',
                                                    dep_var_df=firth_protein_sharing[
                                                        firth_protein_sharing['fdr_sig'] == 1].copy(),
                                                    dep_var_col='n_shared')
# alpha vs SemSim
firth_semsim_vs_prevalent_alpha = linear(indep_var_df=prior_proteomics_sig, indep_var_col='alpha_hat',
                                         dep_var_df=firth_semsim, dep_var_col='bma_sem_sim')
# store results
prior_proteomics_results = {'firth_sig_protein_sharing': firth_sig_protein_sharing_vs_prevalent_alpha,
                            'firth_n_protein_sharing': firth_n_protein_sharing_vs_prevalent_alpha,
                            'firth_semsim': firth_semsim_vs_prevalent_alpha}

""" Prevalent/prior proteomics vs prior-disease-adjusted Firth """
# SemSim
firth_prior_adjusted_semsim = pd.read_csv(
    '~/data/internal/proteomics/pathway_based_similarities/firth_disease_resnik_bma_prior_disease.csv')

# alpha vs SemSim
firth_prior_adjusted_semsim_vs_prevalent_alpha = linear(indep_var_df=prior_proteomics_sig, indep_var_col='alpha_hat',
                                                        dep_var_df=firth_prior_adjusted_semsim,
                                                        dep_var_col='bma_sem_sim')
# store results
prior_proteomics_vs_prior_adjusted_results = {
    'firth_prior_adjusted_semsim': firth_prior_adjusted_semsim_vs_prevalent_alpha}

""" Newly gained co-occurrence in proteomics """
# newly gained co-occurrence
full_proteomics = pd.read_csv('~/data/internal/proteomics/proteomics_pan_ukbb_eur_qu10_50_cooccurrences_alphamle.csv')
full_proteomics = significance_labelling(full_proteomics, pvalue_col='alpha_hat_pval', pair_mirroring=True,
                                         pair_columns=['disease1', 'disease2'])
full_proteomics_sig = full_proteomics[full_proteomics['fdr_sig'] == 1].copy()
full_proteomics_sig = full_proteomics_sig.loc[~full_proteomics_sig.index.duplicated(keep='first'), :]

# to get newly_gained, could use how='left_anti' in pandas 3+
full_proteomics_sig['disease_pair'] = full_proteomics_sig.apply(
    lambda row: '-'.join([row['disease1'], row['disease2']]),
    axis=1)
prior_proteomics_disease_pairs = prior_proteomics_sig.apply(lambda row: '-'.join([row['disease1'], row['disease2']]),
                                                            axis=1).to_list()
newly_gained_proteomics_sig = full_proteomics_sig[
    ~full_proteomics_sig['disease_pair'].isin(prior_proteomics_disease_pairs)].copy()

# protein sharing
cox_protein_sharing = pd.read_csv(
    '~/data/internal/proteomics/shared_assoc_proteins/cox_protein_sharing_fisher_cossims.csv')
cox_protein_sharing = (cox_protein_sharing.merge(
    disease_info[['disease_field', 'code_chapter']]
    .rename(columns={'code_chapter': 'disease1'}))
                       .drop(columns='disease1').rename(columns={'disease_field': 'disease1'})
                       .merge(disease_info[['disease_field', 'code_chapter']]
                              .rename(columns={'code_chapter': 'disease2'}))
                       .drop(columns='disease2').rename(columns={'disease_field': 'disease2'}))

# SemSim
cox_semsim = pd.read_csv(
    '~/data/internal/proteomics/pathway_based_similarities/cox_disease_resnik_bma.csv')

# new alpha vs Cox sig-sharing-or-not
cox_sig_protein_sharing_vs_new_alpha = logistic(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                                dep_var_df=cox_protein_sharing, dep_var_col='fdr_sig')
# new alpha vs Cox number of shared entities
cox_n_protein_sharing_vs_new_alpha = linear(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                            dep_var_df=cox_protein_sharing[cox_protein_sharing['fdr_sig'] == 1].copy(),
                                            dep_var_col='n_shared')
# new alpha vs Cox SemSim
cox_semsim_vs_new_alpha = linear(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                 dep_var_df=cox_semsim, dep_var_col='bma_sem_sim')

# new alpha vs Firth sig-sharing-or-not
firth_sig_protein_sharing_vs_new_alpha = logistic(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                                  dep_var_df=firth_protein_sharing, dep_var_col='fdr_sig')
# new alpha vs Firth number of shared entities
firth_n_protein_sharing_vs_new_alpha = linear(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                              dep_var_df=firth_protein_sharing[
                                                  firth_protein_sharing['fdr_sig'] == 1].copy(),
                                              dep_var_col='n_shared')
# new alpha vs Firth SemSim
firth_semsim_vs_new_alpha = linear(indep_var_df=newly_gained_proteomics_sig, indep_var_col='alpha_hat',
                                   dep_var_df=firth_semsim, dep_var_col='bma_sem_sim')

# store results
newly_gained_proteomics_results = {'cox_sig_protein_sharing': cox_sig_protein_sharing_vs_new_alpha,
                                   'cox_n_protein_sharing': cox_n_protein_sharing_vs_new_alpha,
                                   'cox_semsim': cox_semsim_vs_new_alpha,
                                   'firth_sig_protein_sharing': firth_sig_protein_sharing_vs_new_alpha,
                                   'firth_n_protein_sharing': firth_n_protein_sharing_vs_new_alpha,
                                   'firth_semsim': firth_semsim_vs_new_alpha}

""" Collate results """
collated_results = []
for results_set, results_set_dict in {'genomics': genomics_results, 'prior_proteomics': prior_proteomics_results,
                                      'prior_proteomics_vs_prior_adjusted': prior_proteomics_vs_prior_adjusted_results,
                                      'newly_gained_proteomics': newly_gained_proteomics_results}.items():
    results_dict_list = [{'cooccurrence_omic': results_set, 'test': indiv_result, 'algo': indiv_result_dict['algo'],
                          'coeff': indiv_result_dict['coeff'], 'pval': indiv_result_dict['pvalue']}
                         for indiv_result, indiv_result_dict in results_set_dict.items()]
    collated_results += results_dict_list
collated_results_df = pd.DataFrame(collated_results)
collated_results_df.to_csv('~/data/internal/multimorbidity/cooccurrence_vs_commonality.csv', index=False)
