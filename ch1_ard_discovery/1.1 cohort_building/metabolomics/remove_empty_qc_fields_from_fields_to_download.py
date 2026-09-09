"""
Remove empty metabolomics QC fields. Empty data fields cannot be located by Table Exporter so require removal.
"""
from pandas import read_csv

fields_to_download = read_csv('ukbiobank/dnanexus_data_download/metabolomics/nmr_metabolomics_fields_to_download.txt',
                              header=None)
fields_to_download.columns = ['field']

# returned as fields that could not be found in job-GfK5jXQJYKgfg7Fp53V9VPJ2 in Project1 v2 (project-Gf0zZQQJYKgVf2g0Z2yjf58j)
# these fields could not be found as there are no data in them
empty_qc_fields = ['p23700_i0', 'p23701_i0', 'p23702_i0', 'p23703_i0', 'p23706_i0', 'p23707_i0', 'p23708_i0',
                   'p23711_i0', 'p23712_i0', 'p23713_i0', 'p23714_i0', 'p23715_i0', 'p23716_i0', 'p23718_i0',
                   'p23719_i0', 'p23720_i0', 'p23722_i0', 'p23723_i0', 'p23724_i0', 'p23725_i0', 'p23726_i0',
                   'p23727_i0', 'p23728_i0', 'p23729_i0', 'p23730_i0', 'p23731_i0', 'p23732_i0', 'p23733_i0',
                   'p23739_i0', 'p23740_i0', 'p23741_i0', 'p23771_i0', 'p23777_i0', 'p23778_i0', 'p23781_i0',
                   'p23788_i0', 'p23795_i0', 'p23802_i0', 'p23803_i0', 'p23809_i0', 'p23810_i0', 'p23811_i0',
                   'p23816_i0', 'p23817_i0', 'p23819_i0', 'p23820_i0', 'p23823_i0', 'p23824_i0', 'p23825_i0',
                   'p23826_i0', 'p23827_i0', 'p23828_i0', 'p23830_i0', 'p23831_i0', 'p23832_i0', 'p23837_i0',
                   'p23838_i0', 'p23839_i0', 'p23844_i0', 'p23845_i0', 'p23851_i0', 'p23858_i0', 'p23859_i0',
                   'p23865_i0', 'p23866_i0', 'p23872_i0', 'p23873_i0', 'p23874_i0', 'p23875_i0', 'p23877_i0',
                   'p23899_i0', 'p23905_i0', 'p23906_i0', 'p23909_i0', 'p23912_i0', 'p23914_i0', 'p23919_i0',
                   'p23944_i0', 'p23947_i0']

fields_to_download_no_empty_qc_fields = fields_to_download[~fields_to_download['field'].isin(empty_qc_fields)]
fields_to_download_no_empty_qc_fields.to_csv('ukbiobank/dnanexus_data_download/metabolomics/nmr_metabolomics_fields_to_download_no_empty_qc_fields.txt',
                                             header=False, index=False, lineterminator='\n')