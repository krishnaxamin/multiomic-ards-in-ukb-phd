""" Compare age-of-onsets for certain ARDs to determine which relative 'order'. """
import pandas as pd
import scipy

g20 = 'p131022'
m16 = 'p131870'
h25 = 'p131164'
h26 = 'p131166'
j84 = 'p131528'
m17 = 'p131872'

disease_info = pd.read_csv(
    f"~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv")


def convert_onset_rate_into_onset_counts(time_series_data: pd.DataFrame):
    onset_counts = []
    for data_row in time_series_data.to_dict(orient='records'):
        if data_row['disease_onset_rate'] == 0:
            continue
        onset_counts += [data_row['age']] * int(data_row['num_diagnosed'])
    return onset_counts


""" Import time series """
# genomics
g20_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{g20}_everyone_time_series.csv")
m16_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{m16}_everyone_time_series.csv")
m17_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{m17}_everyone_time_series.csv")
h25_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{h25}_everyone_time_series.csv")
h26_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{h26}_everyone_time_series.csv")
j84_genomics = pd.read_csv(
    f"~/data/internal/time_series/genomics_pan_ukbb_eur_time_series/common_unisex_everyone/{j84}_everyone_time_series.csv")

# proteomics
g20_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{g20}_everyone_time_series.csv")
m16_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{m16}_everyone_time_series.csv")
m17_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{m17}_everyone_time_series.csv")
h25_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{h25}_everyone_time_series.csv")
h26_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{h26}_everyone_time_series.csv")
j84_proteomics = pd.read_csv(
    f"~/data/internal/time_series/proteomics_pan_ukbb_eur_time_series/common_unisex_everyone/{j84}_everyone_time_series.csv")

""" Genomics age-of-onset comparisons """
print(f"G20 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='greater').pvalue}")
print(f"G20 vs H25 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(h25_genomics), alternative='greater').pvalue}")
print(f"H25 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='greater').pvalue}")
print(f"G20 vs H26 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(h26_genomics), alternative='greater').pvalue}")
print(f"H26 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='greater').pvalue}")
print(f"H26 vs M17 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_genomics), convert_onset_rate_into_onset_counts(m17_genomics), alternative='greater').pvalue}")
print(f"M16 vs M17 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_genomics), convert_onset_rate_into_onset_counts(m17_genomics), alternative='greater').pvalue}")
print(f"G20 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='greater').pvalue}")
print(f"M16 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='greater').pvalue}")
print(f"H25 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h25_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='greater').pvalue}")

print(f"G20 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='less').pvalue}")
print(f"G20 vs H25 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(h25_genomics), alternative='less').pvalue}")
print(f"H25 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='less').pvalue}")
print(f"G20 vs H26 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(h26_genomics), alternative='less').pvalue}")
print(f"H26 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_genomics), convert_onset_rate_into_onset_counts(m16_genomics), alternative='less').pvalue}")
print(f"H26 vs M17 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_genomics), convert_onset_rate_into_onset_counts(m17_genomics), alternative='less').pvalue}")
print(f"M16 vs M17 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_genomics), convert_onset_rate_into_onset_counts(m17_genomics), alternative='less').pvalue}")
print(f"G20 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='less').pvalue}")
print(f"M16 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='less').pvalue}")
print(f"H25 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h25_genomics), convert_onset_rate_into_onset_counts(j84_genomics), alternative='less').pvalue}")

""" Proteomics age-of-onset comparisons """
print(f"G20 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='greater').pvalue}")
print(f"G20 vs H25 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(h25_proteomics), alternative='greater').pvalue}")
print(f"H25 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='greater').pvalue}")
print(f"G20 vs H26 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(h26_proteomics), alternative='greater').pvalue}")
print(f"H26 vs M16 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='greater').pvalue}")
print(f"H26 vs M17 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_proteomics), convert_onset_rate_into_onset_counts(m17_proteomics), alternative='greater').pvalue}")
print(f"M16 vs M17 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_proteomics), convert_onset_rate_into_onset_counts(m17_proteomics), alternative='greater').pvalue}")
print(f"G20 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='greater').pvalue}")
print(f"M16 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='greater').pvalue}")
print(f"H25 vs J84 - greater: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h25_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='greater').pvalue}")

print(f"G20 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='less').pvalue}")
print(f"G20 vs H25 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(h25_proteomics), alternative='less').pvalue}")
print(f"H25 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='less').pvalue}")
print(f"G20 vs H26 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(h26_proteomics), alternative='less').pvalue}")
print(f"H26 vs M16 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_proteomics), convert_onset_rate_into_onset_counts(m16_proteomics), alternative='less').pvalue}")
print(f"H26 vs M17 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h26_proteomics), convert_onset_rate_into_onset_counts(m17_proteomics), alternative='less').pvalue}")
print(f"M16 vs M17 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_proteomics), convert_onset_rate_into_onset_counts(m17_proteomics), alternative='less').pvalue}")
print(f"G20 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(g20_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='less').pvalue}")
print(f"M16 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(m16_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='less').pvalue}")
print(f"H25 vs J84 - less: {scipy.stats.mannwhitneyu(convert_onset_rate_into_onset_counts(h25_proteomics), convert_onset_rate_into_onset_counts(j84_proteomics), alternative='less').pvalue}")
