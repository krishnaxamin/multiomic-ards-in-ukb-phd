"""
Get which associated variants are also pQTLs as determined by the UKB-PPP (Sun et al, 2023).
"""
from pandas import read_csv, concat, merge
from liftover import get_lifter
from pyprind import ProgBar
from selenium import webdriver
from selenium.webdriver.common.by import By
from time import sleep

import sys
import requests
import os

server = 'https://rest.ensembl.org'

""" Import and set up data """

geno = read_csv('~/ch2_genomics/2.5 full_gwas_runs/pan-ukbb-eur_assoc_geno_saige_minimal.csv')
imputed = read_csv('~/ch2_genomics/2.5 full_gwas_runs/pan-ukbb-eur_assoc_imputed_saige_minimal.csv')

geno_assoc = geno[geno['assoc_label'] == 'assoc'].copy()
imputed_assoc = imputed[(imputed['assoc_label'] == 'assoc') & (imputed['info_label'] == 'high_info')].copy()

assoc = concat([geno_assoc[['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'UKB_ref', 'UKB_alt', 'Allele2_is_UKB_alt', 'BETA', 'disease']],
                imputed_assoc[['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2', 'UKB_ref', 'UKB_alt', 'Allele2_is_UKB_alt', 'BETA', 'disease']]])
assoc['CHR'] = assoc['CHR'].astype(str)
assoc.drop_duplicates(['CHR', 'POS', 'UKB_ref', 'UKB_alt', 'disease'], inplace=True)
assoc.reset_index(drop=True, inplace=True)

""" Lift over coordinates from hg19 to hg38 where possible """

# liftover converter from hg19 (current assembly) to target assembly of hg38
converter = get_lifter('hg19', 'hg38')

new_chrs = []
new_positions = []
unmappable = []  # 10 unmappable variants
bar = ProgBar(len(assoc), stream=sys.stdout, title='Lifting over')
for i in range(len(assoc)):
    old_chr = assoc.at[i, 'CHR']
    old_pos = assoc.at[i, 'POS']

    new_coords = converter[old_chr][old_pos]
    if len(new_coords) == 0:

        # use Ensembl to get hg38 provided rsID
        rsid = assoc.at[i, 'MarkerID']
        if 'rs' not in rsid:
            unmappable.append(str(old_chr) + ': ' + str(old_pos))
            new_chrs.append('')
            new_positions.append('')
        else:
            ext = '/vep/human/id/' + str(rsid) + '?'
            try:
                r = requests.get(server + ext, headers={"Content-Type": "application/json"})
                r.raise_for_status()
            except requests.exceptions.HTTPError as err:
                print(err.response.text)
                unmappable.append(str(old_chr) + ': ' + str(old_pos))
                new_chrs.append('')
                new_positions.append('')
                continue

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            variant_info = r.json()[0]

            new_chr = variant_info['seq_region_name']
            new_pos = variant_info['start']
            new_chrs.append(new_chr)
            new_positions.append(new_pos)

    else:
        new_chr = new_coords[0][0].split('r')[1]  # this isn't going to change, is it?
        new_pos = new_coords[0][1]
        new_chrs.append(new_chr)
        new_positions.append(new_pos)

    bar.update()

# add hg38 coordinates
assoc['pos_hg38'] = new_positions
# assoc_hg38 = assoc.copy()
# assoc_hg38['CHR'] = new_chrs
# assoc_hg38['POS'] = new_positions

""" Get 200kb ranges encompassing all assoc variants """

assoc_ordered = assoc.sort_values(by=['CHR', 'pos_hg38', 'disease']).reset_index(drop=True)

# some variants will have failed the LiftOver. these variants don't have hg38 coordinates but will be included in the ranges, so if there's an ID match, mapping is possible
assoc_failed_rsids = [x for x in list(set(assoc_ordered[assoc_ordered['pos_hg38'] == '']['MarkerID'])) if 'rs' in x]

# remove failed variants from the main dataset
assoc_ordered = assoc_ordered[assoc_ordered['pos_hg38'] != ''].copy().reset_index(drop=True)

# ranges specified in form [chr, start, end]
chr_list = [str(x) for x in list(range(1, 23)) + ['X']]
range_list = []
for chromo in chr_list:
    assoc_chr = assoc_ordered[assoc_ordered['CHR'] == chromo].copy().reset_index(drop=True)
    while len(assoc_chr) > 0:
        start = list(assoc_chr['pos_hg38'])[0]
        end = start + 200000
        range_list.append([chromo, start, end])
        assoc_chr = assoc_chr[assoc_chr['pos_hg38'] > end].copy()

""" Scrape UKB-PPP portal """

# set up the driver
driver = webdriver.Chrome()

# navigate to UKB-PPP portal
driver.get('https://metabolomips.org/ukbbpgwas/pgwas.table.php#tabs-in-chrpos-range')

# get relevant input boxes and buttons
chr_box = driver.find_element(By.ID, 'in-chrpos-range-chr')
start_box = driver.find_element(By.ID, 'in-chrpos-range-from')
end_box = driver.find_element(By.ID, 'in-chrpos-range-to')
query_button = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[1]/div[6]/div[3]/button')
download_button = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[3]/div/div[1]/div[1]/div/button[1]')

# list to contain query results
query_results = []

# loop through ranges
bar = ProgBar(len(range_list), stream=sys.stdout, title='Querying UKB-PPP')
for query_range in range_list:

    # clear input boxes
    chr_box.clear()
    start_box.clear()
    end_box.clear()

    # fill in input boxes
    chr_box.send_keys(query_range[0])
    start_box.send_keys(query_range[1])
    end_box.send_keys(query_range[2])

    # query
    query_button.click()

    sleep(1)

    # click 'Download' button
    download_button.click()
    sleep(0.2)

    # identify button to download as CSV
    download_csv_button = driver.find_element(By.XPATH, '/html/body/div[1]/div[2]/div[3]/div/div[1]/div[1]/div/div[2]/div/button[2]')

    # download as CSV
    download_csv_button.click()
    sleep(1)

    # click off the download menu to return to normal view
    driver.find_element(By.XPATH, "//html").click()

    # read in downloaded file
    download_dir = ''
    query_result = read_csv(f"{download_dir}\Plasma proteomic associations with genetics and health in the UK Biobank.csv")
    query_results.append(query_result)

    # remove downloaded file
    os.remove(f"{download_dir}\Plasma proteomic associations with genetics and health in the UK Biobank.csv")

    bar.update()

driver.quit()

# concatenate query results
ukb_pqtls = concat(query_results)
ukb_pqtls['chr'] = ukb_pqtls['chr'].astype(str)
ukb_pqtls['pos'] = ukb_pqtls['pos'].str.replace(',', '')
ukb_pqtls['pos'] = ukb_pqtls['pos'].astype(int)
# some pQTLs may have p-values recorded as less-than-10^-324 as a string. this corrects that.
#  1e-324 may appear as 0.0 because it is the smallest positive normalized floating-point number in many double-precision systems
ukb_pqtls.loc[ukb_pqtls['p'].str.contains('10-324'), 'p'] = 1e-324
ukb_pqtls.to_csv('~/data/external/qtls/ukb_ppp_pqtls.csv', index=False)
