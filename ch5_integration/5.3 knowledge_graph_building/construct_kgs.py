"""
Pipeline to build tissue-specific cellular networks using known data and association-determined links.
"""
from pyprind import ProgBar
from liftover import get_lifter
from typing import Union, List, Dict, Tuple, Optional

import networkx as nx
import sys
import os
import requests
import pandas as pd
import numpy as np
import re
import logging
import argparse

########################################################################################################################
#
# CONFIG AREA
#
########################################################################################################################

LOG_PATH = ''

# digest Apocrita job name to give a log file
parser = argparse.ArgumentParser(description='Reading in Apocrita job name.')
parser.add_argument('--jobname', type=str, required=True, help='$JOB_NAME')
args = parser.parse_args()  # parse arguments
job_name = args.jobname
LOGGING_PATH = f"{LOG_PATH}/" + job_name + '.log'

# In a module, __name__ = module's name in the Python package namespace
# __name__ is in the 'name' argument. logger with name foo is the parent of a logger with name foo.bar etc.
logger = logging.getLogger(job_name)
logger.setLevel(logging.INFO)
logger.propagate = False

logger.handlers.clear()  # clear existing handlers, including streaming to console

file_handler = logging.FileHandler(LOGGING_PATH, mode='w')
file_handler.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s - %(name)s. %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)


########################################################################################################################
#
# EXTERNAL CUSTOM FUNCTIONS
#
########################################################################################################################


def current_ensembl_version() -> str:
    """
    Gets current Ensembl version.
    
    :returns: Current Ensembl version, e.g. 115
    :rtype: str
    """
    r = requests.get('https://rest.ensembl.org/info/software', headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
    current_version = str(r.json()['release'])
    return current_version


def generate_ensembl_gene_set(current_version: str) -> pd.DataFrame:
    """
    Pulls the current Ensembl gene set: every gene, its coordinates and its ENSG ID in Ensembl.

    :param current_version: Current ensembl version, e.g. 115
    :type current_version: str
    :returns: Current Ensembl gene set
    :rtype: pd.DataFrame
    """
    ensembl_gene_set = pd.read_csv(
        'https://ftp.ensembl.org/pub/release-' + current_version + '/gtf/homo_sapiens/Homo_sapiens.GRCh38.' + current_version + '.gtf.gz',
        sep='\t', skiprows=5, header=None)
    ensembl_gene_set_only_genes = ensembl_gene_set[ensembl_gene_set[2] == 'gene']
    ensembl_gene_set_only_genes.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                                           'attribute']

    # process attribute column to generate multiple columns for each attribute listed
    regex_pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    attributes = ensembl_gene_set_only_genes['attribute'].map(
        lambda s: dict(regex_pattern.findall(str(s).strip().strip('"').replace('""', '"')))
    ).apply(pd.Series)
    ensembl_gene_set_only_genes = pd.concat([ensembl_gene_set_only_genes, attributes], axis=1)

    # export
    return ensembl_gene_set_only_genes.drop(columns=['score', 'frame', 'attribute'])


def ensg_to_xrefs(ensembl_id, xref: str) -> List[str]:
    """
    Given an Ensembl ID, returns a list of matched IDs from other DBs for a given Ensembl ID. Can return an empty list.
    :param ensembl_id: An Ensembl ID, e.g. ENSG00000123456
    :type ensembl_id: str
    :param xref:
    :type xref: str
    :returns: List of matched IDs for a given Ensembl ID
    :rtype: List[str]
    """

    # ensembl_id = input_tuple[0]

    # ensp_id = 'ENSP00000378058'

    # required to limit rate of API calls
    # time.sleep(0.05)

    server = 'https://rest.ensembl.org'
    ext = '/xrefs/id/' + ensembl_id + '?species=human'

    r = requests.get(server + ext, headers={'Content-Type': 'application/json'})

    if not r.ok:
        try:
            r.raise_for_status()
            # sys.exit()
        except requests.exceptions.HTTPError as e:
            print(f"{ensembl_id} does not have Ensembl information.")
            return {}

    xref_info = r.json()

    # matches to UniProt
    if xref == 'uniprot':
        xref_matches = [xref['primary_id'] for xref in xref_info if 'Uniprot' in xref['dbname']]
    elif xref == 'hugo':
        xref_matches = [xref['primary_id'] for xref in xref_info if xref['dbname'] == 'HGNC']
    elif xref == 'ncbi':
        xref_matches = [xref['primary_id'] for xref in xref_info if xref['dbname'] == 'EntrezGene']
    else:
        logger.info(f"{xref} is not covered in the API querying of ENSG to {xref}.")
        return []

    return xref_matches


def propagation_up_ontology_hierarchy(id_ancestor_links: pd.DataFrame, annotations_to_propagate: pd.DataFrame,
                                      annotation_col_id: str, relation_class: Optional[str] = None,
                                      relation_class_col_id: Optional[str] = None) -> pd.DataFrame:
    """
    Propagate annotations up a provided ontology tree, such that entities annotated with a child term are also annotated
     with the parent term.
    :param id_ancestor_links: DataFrame describing child-ancestor (child-parent, child-grandparent, etc.) links.
    Must have columns 'id' and 'ancestor'
    :type id_ancestor_links: pd.DataFrame
    :param annotations_to_propagate: DataFrame of annotations to propagate
    :type annotations_to_propagate: pd.DataFrame
    :param annotation_col_id: column in annotations_to_propagate containing the annotation terms
    :type annotation_col_id: str
    :param relation_class: type of relationship between child and ancestor
    :type relation_class: Optional[str]
    :param relation_class_col_id: column in annotations_to_propagate containing relation_class
    :type relation_class_col_id: Optional[str]
    :return: original annotations with annotations fully propagated up the ontology tree
    :rtype: pd.DataFrame
    """
    # argument check
    if relation_class is not None:
        assert relation_class_col_id is not None, 'To specify a relation_class, the col_id must also be specified.'
    if relation_class is not None:
        id_ancestor_links = id_ancestor_links[id_ancestor_links[relation_class_col_id] == relation_class].copy()
    else:
        id_ancestor_links = id_ancestor_links.copy()

    annotations_to_propagate_list = list(annotations_to_propagate[annotation_col_id].unique())
    filtered_links = id_ancestor_links[id_ancestor_links['id'].isin(annotations_to_propagate_list)].copy()

    # generate mapping: {id: [list of ancestors, including id itself]}
    id_to_ancestor_map = (filtered_links
                          .groupby('id')['ancestor']
                          .unique()
                          .apply(list)
                          .to_dict())

    annotations_to_propagate_with_ancestors = annotations_to_propagate.copy()
    annotations_to_propagate_with_ancestors['ancestors'] = annotations_to_propagate_with_ancestors[
        annotation_col_id].apply(
        lambda x: id_to_ancestor_map[x])

    annotations_propagated = (annotations_to_propagate_with_ancestors.explode('ancestors', ignore_index=True)
                              .drop(columns=annotation_col_id)
                              .rename(columns={'ancestors': annotation_col_id})
                              .drop_duplicates(ignore_index=True)
                              .merge(annotations_to_propagate.assign(original_annotation=1),
                                     on=list(annotations_to_propagate.columns),
                                     how='outer')
                              .fillna(0)
                              .astype({'original_annotation': 'int8'}))

    return annotations_propagated


########################################################################################################################
#
# READ IN DATA
#
########################################################################################################################


def read_in_data() -> Dict[str, Union[List[str], pd.DataFrame, str]]:
    """
    Loads in data to be used in graph construction from various sources.
    :return: Dictionary containing data to be used in graph construction.
    :rtype: Dict[str, Union[List[str], pd.DataFrame, str]]
    """
    logger.info('Starting to read in data.')
    # read in GWAS hits
    hits = pd.read_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_lifestyles.csv', dtype={'CHR': str})
    hits = hits[(hits['assoc_label'] == 'assoc') & (hits['info_label'] == 'high_info')].reset_index(drop=True).copy()
    hits = hits[['CHR', 'POS', 'MarkerID', 'UKB_ref', 'UKB_alt', 'BETA', 'Allele2_is_UKB_alt', 'disease']].rename(
        columns={'BETA': 'effect'}).copy()
    logger.info('GWAS hits read in successfully.')

    # get file names for eQTLs from each tissue
    eqtl_files = [x.name for x in os.scandir('~/data/external/qtls/GTEx_Analysis_v8_eQTL') if
                  x.is_file() and 'signif' in x.name]
    gtex_tissues = sorted([x.split('.')[0] for x in eqtl_files])

    eqtl_paths = ['~/data/external/qtls/GTEx_Analysis_v8_eQTL/' + x for x in eqtl_files]
    logger.info('Paths to GTEx eQTL data determined successfully.')

    # read in disease info
    disease_info = pd.read_csv('~/data/internal/phenotype_coding/pan_ukbb_eur_qu10_50_diseases_icd10-info_num-diagnosed-info.csv')
    disease_info['code_chapter'] = disease_info['icd10_three_letter'] + ' (' + disease_info['icd10_chapter'].astype(
        str) + ')'
    logger.info('Disease information read in successfully.')

    # read in pQTLs
    pqtls = pd.read_csv('~/data/external/qtls/ukb_ppp_pqtls.csv')[
        ['chr', 'pos', 'rsid', 'Target.1', 'beta']].rename(columns={'beta': 'effect'})
    logger.info('pQTL data from the UKB-PPP read in succesfully.')

    # read in gene-disease MAGMA
    magma = pd.read_csv('~/data/internal/genomics/magma/lifestyles/gene_level_results_assoc.csv')
    magma = magma[magma['bonf_sig'] == 1].copy()
    magma = magma.merge(
        disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(
        columns='disease').rename(columns={'icd10_three_letter': 'disease'})[['GENE', 'disease']].copy()

    # read in Firth proteomics
    firth = pd.read_csv('~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_firth.csv')
    firth = firth.merge(
        disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(
        columns='disease').rename(columns={'icd10_three_letter': 'disease', 'beta': 'effect'})
    firth = firth[firth['fdr_sig'] == 1].copy()
    logger.info('Firth proteomics data read in successfully.')

    # read in Cox proteomics
    cox = pd.read_csv('~/data/internal/proteomics/pan-ukbb-eur_assoc_proteomics_cox.csv')
    cox = cox.merge(
        disease_info[['disease_field', 'icd10_three_letter']].rename(columns={'disease_field': 'disease'})).drop(
        columns='disease').rename(columns={'icd10_three_letter': 'disease', 'beta': 'effect'})
    cox = cox[cox['fdr_sig'] == 1].copy()
    logger.info('Cox proteomics data read in successfully.')

    # read in stringdb
    stringdb = pd.read_csv('~/data/external/stringdb/9606.protein.links.full.v12.0.txt.gz', sep='\s+')
    # isolate high-confidence stringdb
    stringdb_coexp = stringdb[stringdb['coexpression'] >= 700][['protein1', 'protein2', 'coexpression']].copy()
    stringdb_int = stringdb[stringdb['experiments'] >= 700][['protein1', 'protein2', 'experiments']].copy()
    logger.info('StringDB high-confidence co-expression and physical interaction links read in successfully.')

    # gene ontology
    go_annotations = pd.read_csv('~/data/external/gene_ontology/goa_human.gaf.gz', sep='\t',
                                 comment='!', header=None, low_memory=False)
    go_annotations.columns = ['id_db', 'id', 'symbol', 'relation', 'go_id', 'reference', 'evidence', 'additional_id',
                              'aspect', 'name', 'synonym', 'type', 'taxon', 'annotation_date', 'annotation_by',
                              'annotation_extension', 'gene_product_id']
    go_cc = go_annotations[(go_annotations['type'] == 'protein') &
                           (~go_annotations['evidence'].isin(['IEA', 'NAS', 'ND', 'NR'])) &
                           (~go_annotations['relation'].str.contains('NOT', na=False)) &
                           (go_annotations['aspect'] == 'C')][
        ['id', 'go_id']].drop_duplicates(ignore_index=True)
    # propagate annotations up hierarchy
    go_cc = propagation_up_ontology_hierarchy(id_ancestor_links=pd.read_csv(
        GROUP_APOCRITA_DIR + '/gene_ontology/go_id_ancestor_relations.csv'),
        annotations_to_propagate=go_cc,
        annotation_col_id='go_id')
    # join multiple GO-CC for one gene into one string
    go_cc = (
        go_cc.groupby('id')['go_id']
        .agg(lambda x: '|'.join(map(str, sorted(set(x)))))
        .reset_index()
    ).drop_duplicates(ignore_index=True)
    logger.info(
        'Gene Ontology Cellular Component annotations read in, with annotations propagated up the hierarcy, successfully.')

    # HPA data
    hpa_exp = pd.read_csv('~/data/external/human_protein_atlas/rna_tissue_detail_gtex.tsv', sep='\t')
    hpa_tissue_specificity = pd.read_csv('~/data/external/human_protein_atlas/tissue_specificity_with_tissue_id.csv')
    hpa_secretome = pd.read_csv('~/data/external/human_protein_atlas/sa_location_Secreted.tsv', sep='\t')
    hpa_secretome = (hpa_secretome[['Gene', 'Secretome location', 'Secretome function']]
                     .set_axis(['Gene name', 'secretome_location', 'secreted_function'], axis=1)
                     .merge(hpa_exp[['Gene', 'Gene name']].drop_duplicates())
                     .drop(columns='Gene name')).drop_duplicates(ignore_index=True)
    logger.info(
        'HPA/GTEx data (tissue-specific expression, specificity categories, secretome function and localisation) read in successfully.')

    logger.info('Data read complete successfully.')

    # construct data_dict
    data_dict = {'hits': hits,
                 'eqtl_paths': eqtl_paths,
                 'gtex_tissues': gtex_tissues,
                 'disease_info': disease_info,
                 'pqtls': pqtls,
                 'magma': magma,
                 'firth': firth,
                 'cox': cox,
                 'stringdb_coexp': stringdb_coexp,
                 'stringdb_int': stringdb_int,
                 'go_cc': go_cc,
                 'hpa_exp': hpa_exp,
                 'hpa_tissue_specificity': hpa_tissue_specificity,
                 'hpa_secretome': hpa_secretome,
                 # 'current_ensembl_version': current_ensembl_version(),
                 'current_ensembl_version': '115'
                 }

    return data_dict


########################################################################################################################
#
# PROCESS DATA
#
########################################################################################################################


def map_hg19_to_hg38(hg19_data: pd.DataFrame, chr_col: str = 'CHR', pos_col: str = 'POS') -> pd.DataFrame:
    """
    Maps a set of single-point genetic features (e.g., variants) from hg19 to hg38. Retains hg19 positional data, and
     prints to log any variants that were not mapped successfully.
    :param hg19_data: DataFrame of variant GWAS hits on hg19
    :type hg19_data: pd.DataFrame
    :param chr_col: column containing variant chromosomes in hg19_data
    :type chr_col: str
    :param pos_col: column containing variant positions in hg19_data
    :type pos_col: str
    :return: DataFrame of variant GWAS hits on hg38, retaining hg19 coordinates
    :rtype: pd.DataFrame
    """
    logger.info('Mapping hg19 variants to hg38 variants using liftover and the Ensembl API.')

    server = 'https://rest.ensembl.org'

    hg19_data_working = hg19_data.copy()

    # liftover converter from hg19 (current assembly) to target assembly of hg38
    converter = get_lifter('hg19', 'hg38')

    new_chrs = []
    new_positions = []
    unmappable = []  # 10 unmappable variants
    # bar = ProgBar(len(hg19_data), stream=sys.stdout, title='Lifting over')
    for i in range(len(hg19_data)):
        old_chr = hg19_data.at[i, chr_col]
        old_pos = hg19_data.at[i, pos_col]

        new_coords = converter[old_chr][old_pos]
        if len(new_coords) == 0:

            # use Ensembl to get hg38 provided rsID
            rsid = hg19_data.at[i, 'MarkerID']
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
                new_chrs.append(str(new_chr))
                new_positions.append(int(new_pos))

        else:
            new_chr = new_coords[0][0].split('r')[1]  # this isn't going to change, is it?
            new_pos = new_coords[0][1]
            new_chrs.append(str(new_chr))
            new_positions.append(int(new_pos))

        # bar.update()

    hg19_data_working['pos_hg19'] = hg19_data_working[pos_col]
    hg19_data_working[pos_col] = new_positions

    hg19_data_working['unique_id'] = (hg19_data_working[chr_col].astype(str) + '_' +
                                      hg19_data_working[pos_col].astype(str) + '_' +
                                      hg19_data_working['UKB_ref'] + '_' + hg19_data_working['UKB_alt'])

    # remove unmapped variants while recording which ones were unmapped
    hg19_data_working_all_mapped = hg19_data_working[hg19_data_working['POS'] != ''].copy()
    hg19_data_working_all_mapped['POS'] = hg19_data_working_all_mapped['POS'].astype(int)

    # record which variants were unmapped
    unmapped_variants = hg19_data_working[hg19_data_working['POS'] == '']['MarkerID'].to_list()
    unmapped_variants_reporter = f"{len(unmapped_variants)} variants were not mappable:\n"
    for unmapped_variant in unmapped_variants:
        unmapped_variants_reporter += f" - {unmapped_variant}\n"
    logger.warning(unmapped_variants_reporter)

    logger.info(
        f"Mapping variants from hg19 to hg38 completed successfully, with {len(unmapped_variants)} variants unmapped.")

    return hg19_data_working_all_mapped  # transformed to have hg38 data in pos_col, added uniqueID based on hg38 coords


def map_variants_to_eqtls(variant_data: pd.DataFrame, eqtl_path: str, eqtl_tissue: str,
                          variant_chr_col: str = 'CHR', variant_pos_col: str = 'POS') -> pd.DataFrame:
    """
    Maps variants on hg38 to GTEx v8 eQTLs, given a path to a specific eQTL dataset.
    :param variant_data: DataFrame containing variant positional data
    :type variant_data: pd.DataFrame
    :param eqtl_path: path to eQTL dataset. eQTL dataset must have cols 'variant_id', 'gene_id', 'slope' and 'pval_beta'
    :type eqtl_path: str
    :param eqtl_tissue: tissue to which these eQTLs are specific
    :type eqtl_tissue: str
    :param variant_chr_col: column containing variant chromosomes in variant_data
    :type variant_chr_col: str
    :param variant_pos_col: column containing variant positions in variant_data
    :type variant_pos_col: str
    :return: DataFrame containing variant IDs that are eQTLs (as CHR_POS_REF_ALT) ('unique_id'), their eGenes (as ENSG)
     ('gene_id'), and that relationship's effect size ('effect')
    """

    logger.info('Mapping variants to tissue-specific GTEx eQTLs')

    # for i, eqtl_path in enumerate(paths_to_eqtl_data):
    eqtls = pd.read_csv(eqtl_path, sep='\t')[['variant_id', 'gene_id', 'slope', 'pval_beta']]
    eqtls[['CHR', 'POS', 'GTEx_ref', 'GTEx_alt', 'build']] = eqtls['variant_id'].str.split('_', expand=True)
    eqtls['POS'] = eqtls['POS'].astype(int)
    eqtls['CHR'] = [x.replace('chr', '') for x in eqtls['CHR'].to_list()]

    hit_eqtls = (eqtls.rename(columns={'slope': 'effect'})
                 .merge(variant_data.drop(columns='effect')
                        .rename(columns={variant_chr_col: 'CHR', variant_pos_col: 'POS'}),
                        on=['CHR', 'POS']))

    hit_eqtls['gene_id'] = [x.split('.')[0] for x in hit_eqtls['gene_id'].to_list()]
    hit_eqtls = hit_eqtls[['unique_id', 'gene_id', 'effect']].drop_duplicates(ignore_index=True)

    # export
    # os.makedirs(f"{SCRATCH_APOCRITA_DIR}qtl_mappings/", exist_ok=True)
    hit_eqtls.to_csv(f"~/data/internal/knowledge_graph/qtl_mappings/assoc_vars_eqtls_{eqtl_tissue}.csv",
                     index=False)

    logger.info(f"Mapping variants to eQTLs complete successfully and exported to "
                f"'~/data/internal/knowledge_graph/qtl_mappings/assoc_vars_eqtls_{eqtl_tissue}.csv'")

    return hit_eqtls


def expand_row(row: pd.Series, col1: str, col2: str, col3: str) -> pd.DataFrame:
    """
    Deal with peculiarities in mapping ENSGs-UniProts-Olink protien fields, where a single protein field might match to
     multiple genes, or vice versa.
    :param row: One row of the current mapping DataFrame
    :type row: pd.Series
    :param col1: name of the first column in 'row' to look for multiplicity
    :type col1: str
    :param col2: name of the second column in 'row' to look for multiplicity
    :type col2: str
    :param col3: name of the column to keep invariant and spread out over expanded multiplicity from the other columns
    :type col3: str
    :return: mapping DataFrame with one ENSG, one UniProt and one protein field
    :rtype: pd.DataFrame
    """
    c1, c2, c3 = row[col1], row[col2], row[col3]

    s1, s2 = c1.split("_"), c2.split("_")

    if "_" in c1 and "_" in c2 and len(s1) == len(s2):
        # Case 1: same number of underscores → pair positionally
        return pd.DataFrame({col1: s1, col2: s2, col3: [c3] * len(s1)})

    elif "_" not in c1 and "_" not in c2:
        # Case 2: no underscores → keep as is
        return pd.DataFrame({col1: [c1], col2: [c2], col3: [c3]})

    else:
        # Case 3: only one column has underscores → expand only that one
        if "_" in c1:
            return pd.DataFrame({col1: s1, col2: [c2] * len(s1), col3: [c3] * len(s1)})
        else:
            return pd.DataFrame({col1: [c1] * len(s2), col2: s2, col3: [c3] * len(s2)})


def pull_current_ensembl_gene_set(current_ensembl_version: str) -> pd.DataFrame:
    """
    Pull the current Enesmbl gene set, either from a current local file (provided it is the current version) or from
     Ensembl, should the local version not exist or be outdated.
    :param current_ensembl_version: The current Ensembl version, e.g. 115
    :type current_ensembl_version: str
    :return: The current Ensembl gene set
    :rtype: pd.DataFrame
    """
    logger.info(f"Pulling current Ensembl gene set: version {current_ensembl_version}")

    ensembl_gene_set_path = '~/data/external/ensembl/ensembl_gene_set_hg38_v' + current_ensembl_version + '.csv'
    if os.path.exists(ensembl_gene_set_path):
        ensembl_gene_set = pd.read_csv('~/data/external/ensembl/ensembl_gene_set_hg38_v' + current_ensembl_version + '.csv')
        logger.info('Current Ensembl gene set read in from previous file.')
    else:
        ensembl_gene_set = generate_ensembl_gene_set(current_version=current_ensembl_version)
        ensembl_gene_set.to_csv(ensembl_gene_set_path, index=False)
        logger.info(f"Current Ensembl gene set generated and exported to {ensembl_gene_set_path}.")

    return ensembl_gene_set


def generate_ensg_uniprot_mappings(current_ensembl_version: str, ensembl_gene_set: pd.DataFrame) -> pd.DataFrame:
    """
    Identifies UniProt IDs for all ENSGs in the given gene set. Pulls either from an existing local file or queries
     the Ensembl API should the local version not exist or be outdated.
    :param current_ensembl_version: The current Ensembl version, e.g. 115
    :type current_ensembl_version: str
    :param ensembl_gene_set: The current Ensembl gene set
    :type ensembl_gene_set: pd.DataFrame
    :return: Mapping DataFrame mapping single ENSG IDs to all matching UniProts (one ENSG - many UniProts)
    :rtype: pd.DataFrame
    """
    logger.info('Generating mappings from ENSG IDs to UniProt IDs.')

    # UniProt to ENSG mapping - uniprot, xref, ensg_id
    ensg_uniprot_path = '~/data/external/ensembl/ensembl_genes_v' + current_ensembl_version + '_uniprot.csv'
    if os.path.exists(ensg_uniprot_path):
        uniprot_xref = pd.read_csv(ensg_uniprot_path).drop_duplicates(ignore_index=True)
        logger.info('Mapping already exists.')
    else:
        # ensembl_gene_set = pd.read_csv(ensembl_gene_set_path)
        ensgs = ensembl_gene_set['gene_id'].to_list()

        uniprot_xref = pd.DataFrame()
        # bar = ProgBar(len(ensgs), stream=sys.stdout, title='Getting xrefs for Ensembl genes')
        for ensg in ensgs:
            uniprot_matches = ensg_to_xrefs(ensg, xref='uniprot')
            if len(uniprot_matches) != 0:
                uniprot_xref = (pd.concat(
                    [uniprot_xref,
                     pd.DataFrame({'ensg_id': [ensg] * len(uniprot_matches), 'uniprot': uniprot_matches})])
                                .drop_duplicates(ignore_index=True))

            # bar.update()

        uniprot_xref.to_csv(ensg_uniprot_path, index=False)
        logger.info(f"Mapping successfully generated and exported to {ensg_uniprot_path}.")

    return uniprot_xref


def id_mappings_for_olink_proteins(
        uniprot_xref: pd.DataFrame) -> pd.DataFrame:  # if the path to result doesn't exist, run this, otherwise do not
    """
    Generate mapping DataFrame for Olink protein field-UniProt-ENSG ID.
    :param uniprot_xref: DataFrame mapping ENSG IDs to UniProts
    :type uniprot_xref: pd.DataFrame
    :return: DataFrame mapping Olink protein fields to UniProts and ENSG IDs
    :rtype: pd.DataFrame
    """

    logger.info('Generating mappings between different IDs for Olink proteins specifically.')

    # check for result presence
    export_file_path = '~/data/internal/proteomics/olink_field_uniprot_ensg_genesymbol_mapping.csv'

    # if mappings already exist, then the function quits here essentially
    if os.path.exists(export_file_path):
        uniprot_genesymbol_ensg_expanded = pd.read_csv(export_file_path).drop_duplicates(ignore_index=True)
        logger.info('Mapping already exists.')
        return uniprot_genesymbol_ensg_expanded

    else:  # else, making the mapping files
        # map olink fields and their uniprot to ENSGs and gene symbols

        uniprot_genesymbol = pd.read_csv('~/data/external/proteomics/olink_assay.dat', sep='\t')[
            ['Assay', 'UniProt']].rename(
            columns={'Assay': 'gene_symbol', 'UniProt': 'uniprot'})
        uniprot_genesymbol['protein_field'] = uniprot_genesymbol.gene_symbol.str.lower()
        uniprot_genesymbol['protein_field'] = uniprot_genesymbol.protein_field.replace(
            {'hla-dra': 'hla_dra', 'hla-a': 'hla_a', 'hla-e': 'hla_e', 'ervv-1': 'ervv_1'})
        uniprot_genesymbol_expanded = pd.concat(
            [expand_row(row, col1='gene_symbol', col2='uniprot', col3='protein_field')
             for _, row in uniprot_genesymbol.iterrows()], ignore_index=True)

        uniprot_genesymbol_ensg_expanded = (uniprot_genesymbol_expanded.merge(uniprot_xref)
                                            .drop_duplicates(ignore_index=True))

        # os.makedirs(HOME_APOCRITA_DIR + 'proteomics/', exist_ok=True)
        uniprot_genesymbol_ensg_expanded.to_csv(export_file_path,
                                                index=False)  # uniprot, ensg_id, gene_symbol, protein_field
        logger.info(f"Mapping successfully generated and exported to {export_file_path}.")

    return uniprot_genesymbol_ensg_expanded


def map_variants_to_pqtls(variant_data: pd.DataFrame, pqtl_data: pd.DataFrame, olink_uniprot_ensg_mapping: pd.DataFrame,
                          variant_chr_col: str = 'CHR', variant_pos_col: str = 'POS') -> pd.DataFrame:
    """
    Maps variants on hg38 to UKB-PPP pQTLs.
    :param variant_data: DataFrame containing variant positional data
    :type variant_data: pd.DataFrame
    :param pqtl_data: DataFrame of UKB-PPP pQTLs. Must have columns 'chr', 'pos', 'rsid', and 'Target.1'
    :type pqtl_data: pd.DataFrame
    :param olink_uniprot_ensg_mapping: DataFrame mapping Olink protein fields to UniProt and ENSG IDs
    :type olink_uniprot_ensg_mapping: pd.DataFrame
    :param variant_chr_col: Column containing variant chromosomes in variant_data
    :type variant_chr_col: str
    :param variant_pos_col: Column containing variant positions in variant_data
    :type variant_pos_col: str
    :return: DataFrame containing variant IDs that are pQTLs (as CHR_POS_REF_ALT) ('unique_id'), their pGenes (as ENSG)
     ('ensg_id'), and that relationship's effect size ('effect')
     :rtype: pd.DataFrame
    """
    logger.info('Mapping variants to UKB-PPP pQTLs.')

    # map by position
    mapping_pos = variant_data.drop(columns='effect').merge(pqtl_data, left_on=[variant_chr_col, variant_pos_col],
                                                            right_on=['chr', 'pos'])

    # map by rsID
    mapping_ids = variant_data.drop(columns='effect').merge(pqtl_data, left_on='MarkerID', right_on='rsid')

    # combine mappings
    mapping = pd.concat([mapping_pos, mapping_ids]).drop_duplicates().sort_values(
        by=['disease', variant_chr_col, variant_pos_col]).reset_index(drop=True)

    # adjust column names
    mapping.rename(columns={'chr': 'chr_from_ukbppp', 'pos': 'pos_from_ukbppp', 'MarkerID': 'id_from_gwas',
                            'rsid': 'id_from_ukbppp', 'Target.1': 'uniprot'},
                   inplace=True)

    # add ENSG IDs
    mapping = mapping.merge(olink_uniprot_ensg_mapping[['uniprot', 'ensg_id']], on='uniprot').drop_duplicates()

    # export
    os.makedirs('~/data/internal/knowledge_graph/qtl_mappings/', exist_ok=True)
    mapping.to_csv('~/data/internal/knowledge_graph/qtl_mappings/assoc_vars_ukbppp_pqtls.csv', index=False)

    logger.info(f"Mapping variants to pQTLs complete successfully and exported to "
                f"'~/data/internal/knowledge_graph/qtl_mappings/assoc_vars_ukbppp_pqtls.csv'")

    return mapping[
        ['unique_id', 'ensg_id', 'effect']].drop_duplicates()  # only returns the variant ID and the QTL'ed ENSG ID


def containing_gene(variant_pos: Union[str, int], gene_set_df_by_chromo: pd.DataFrame) -> Tuple[str, int]:
    """
    For a given variant, identify the closest gene(s).
    :param variant_pos: The position of the variant
    :type variant_pos: Union[str, int]
    :param gene_set_df_by_chromo: Ensembl gene set genes from the same chromosome as the variant.
    :type gene_set_df_by_chromo: pd.DataFrame
    :return: Gene ENSG ID, and the distance from the variant to that gene. Distance set to -10 if there is no proximal
     gene. Distance can be 0 if the variant is in a gene. Multiple genes can be returned, with IDs concatenated by '|'.
    :rtype: Tuple[str, int]
    """
    # check if variant is in a gene
    variant_in_genes = gene_set_df_by_chromo[(gene_set_df_by_chromo['start'] <= int(variant_pos)) &
                                             (gene_set_df_by_chromo['end'] >= int(variant_pos))].copy()
    if not variant_in_genes.empty:
        proximal_genes = '|'.join(variant_in_genes['gene_id'].to_list())
        dist = 0
        return proximal_genes, dist
    else:
        return '', -10


def variant_gene_proximity_mapping(variant_data: pd.DataFrame, ensembl_gene_set_working: pd.DataFrame,
                                   variant_chr_col: str = 'CHR', variant_pos_col: str = 'POS') -> pd.DataFrame:
    """
    Find the closest genes to a given set of hg38 variants.
    :param variant_data: DataFrame of variants in hg38 coordinates.
    :type variant_data: pd.DataFrame
    :param ensembl_gene_set_working: Current Ensembl geen set.
    :type ensembl_gene_set_working: pd.DataFrame
    :param variant_chr_col: Column containing variant chromosomes in variant_data
    :type variant_chr_col: str
    :param variant_pos_col: Column containing variant positions in variant_data
    :type variant_pos_col: str
    :return: DataFrame mapping variants in variant_data to closest ENSG IDs, in a 1:1 manner
    :rtype: pd.DataFrame
    """
    logger.info('Identifying genes most proximal to variants.')

    # only concerned with genes containing variants
    ensembl_gene_set_working_func = ensembl_gene_set_working[['seqname', 'start', 'end', 'gene_id']].copy()
    ensembl_gene_set_working_func.sort_values(by=['seqname', 'start'], inplace=True, ignore_index=True)
    # split the gene_set by chromosome
    ensembl_gene_sets_by_chromo = {chromo: chromo_df
                                   for chromo, chromo_df in ensembl_gene_set_working_func.groupby('seqname')}

    variant_data_working = variant_data.copy()
    variant_data_working[['proximal_gene', 'dist_to_proximal_gene']] = variant_data_working.apply(
        lambda x: pd.Series(containing_gene(x[variant_pos_col], ensembl_gene_sets_by_chromo[x[variant_chr_col]])),
        axis=1)

    variant_data_in_genes = variant_data_working[variant_data_working['proximal_gene'] != ''][
        ['unique_id', 'proximal_gene', 'dist_to_proximal_gene']].rename(columns={'proximal_gene': 'gene_id'}).copy()

    variant_data_in_genes['gene_id'] = variant_data_in_genes['gene_id'].str.split('|')
    variant_data_in_genes = variant_data_in_genes.explode('gene_id').reset_index(drop=True).drop_duplicates()

    logger.info('Variant-to-proximal-gene mapping completed successfully.')

    return variant_data_in_genes  # returns only variant ID and ENSG ID for those variants contained within a gene


def make_gwas_edges(variant_data: pd.DataFrame) -> pd.DataFrame:
    """
    Make DataFrame for Variant -> Disease GWAS edges.
    :param variant_data: DataFrame of variant data in hg38 coordinates.
    :type variant_data: pd.DataFrame
    :return: DataFrame mapping variants to disease in a 1:1: manner, with effect size ('effect')
    :rtype: pd.DataFrame
    """
    logger.info('Making Variant -> Disease edges from GWAS data.')

    variant_data_working = variant_data[['unique_id', 'effect', 'Allele2_is_UKB_alt', 'disease']].copy()
    # flip betas if UKB_ref is used as Allele2 in SAIGE to calculate effect size, so that UKB_alt is now Allele2
    variant_data_working.loc[variant_data_working['Allele2_is_UKB_alt'], 'effect'] = -1 * variant_data_working.loc[
        variant_data_working['Allele2_is_UKB_alt'], 'effect']

    logger.info('Variant -> Disease edges created successfully.')

    return variant_data_working.drop(columns='Allele2_is_UKB_alt').drop_duplicates()  # unique_id, effect, disease


def make_proteomics_edges(
        proteomics_data: pd.DataFrame, olink_ensg_mapping: pd.DataFrame, algo: str = 'none') -> pd.DataFrame:
    """
    Make DataFrame for Gene -> Disease proteomics edges.
    :param proteomics_data: DataFrame of proteomics results
    :type proteomics_data: pd.DataFrame
    :param olink_ensg_mapping: DataFrame mapping Olink protein fields to UniProts and ENSG IDs
    :type olink_ensg_mapping: pd.DataFrame
    :param algo: specify whether it's Firth or Cox (or other) results the edges are being made to represent
    :type algo: str
    :return: DataFrame mapping genes to disease in a 1:1 manner, with effect size ('effect')
    """
    logger.info('Making Gene -> Disease edges from proteomics data.')
    # call function for firth and cox data together separately
    proteomics_data_working = proteomics_data[['protein', 'disease', 'effect']].rename(
        columns={'protein': 'protein_field'}).copy()
    proteomics_data_working = (proteomics_data_working
                               .merge(olink_ensg_mapping[['protein_field', 'ensg_id']])
                               .drop(columns='protein_field'))

    # edge case: two protein_fields for one UniProt/ENSG -> two effects for one ENSG -> merge by averaging
    proteomics_data_working_averaged = (
        proteomics_data_working.groupby(['ensg_id', 'disease'], as_index=False)
        .agg({'effect': lambda x: np.mean(x)})
    ).drop_duplicates()

    logger.info(f"Gene -> Disease edges created successfully (algo: {algo}).")

    return proteomics_data_working_averaged  # ensg_id, disease, effect


def make_magma_edges(current_ensembl_version: str, magma_df: pd.DataFrame) -> pd.DataFrame:
    """
    Make DataFrame for Gene -> Disease MAGMA edges.
    :param current_ensembl_version: The current Ensembl version, e.g. 115
    :type current_ensembl_version: str
    :param magma_df: DataFrame of Bonferroni-significant MAGMA-derived gene-disease associations.
    :type magma_df: pd.DataFrame
    :return:
    """
    ensembl_gene_set = pd.read_csv(
        f"~/data/external/ensembl/ensembl_gene_set_hg38_v{current_ensembl_version}.csv")

    ensgs = ensembl_gene_set['gene_id'].to_list()

    # get xrefs
    ncbi_ensembl_mapping = pd.DataFrame()
    bar = ProgBar(len(ensgs), stream=sys.stdout, title='Getting xrefs for Ensembl genes')
    for ensg in ensgs:
        ncbi_matches = ensg_to_xrefs(ensg, xref='ncbi')
        if len(ncbi_matches) != 0:
            ncbi_ensembl_mapping = pd.concat(
                [ncbi_ensembl_mapping, pd.DataFrame({'ensg_id': [ensg] * len(ncbi_matches), 'ncbi': ncbi_matches})])

        bar.update()

    magma_mapped = (magma_df
                    .rename(columns={'GENE': 'ncbi'})
                    .merge(ncbi_ensembl_mapping)
                    .drop(columns='ncbi')
                    .drop_duplicates(ignore_index=True))

    logger.info('Gene -> Disease edges (MAGMA) created successfull.')

    return magma_mapped  # ensg_id, disease


def annotate_genes(ensembl_genes: pd.DataFrame, go_cc_annots: pd.DataFrame, tissue_exp: pd.DataFrame,
                   tissue_specificity: pd.DataFrame, secretome_data: pd.DataFrame,
                   ensg_uniprot_mapping: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Annotate genes with (1) GO-CC annotations (2) tissue-specific expression info (from GTEx/HPA) (3) secretome info
     (from HPA). The gene set is then split based on tissue, and gnes with nTPM <= 0.1 are taken as not expressed in
     that tissue.
    :param ensembl_genes: Current Ensembl gene set
    :type ensembl_genes: pd.DataFrame
    :param go_cc_annots: DataFrame mapping UniProt IDs to GO-CC annotations. One row per UniProt, one UniProt can map
     to many GO-CC.
    :type go_cc_annots: pd.DataFrame
    :param tissue_exp: DataFrame providing nTPM expression data for GTEx tissues for ENSG IDs
    :type tissue_exp: pd.DataFrame
    :param tissue_specificity: DataFrame providing tissue-specific expression labels in GTEx tissues for ENSG IDs
    :type tissue_specificity: pd.DataFrame
    :param secretome_data: DataFrame providing secretome data from HPA for ENSG IDs
    :type secretome_data: pd.DataFrame
    :param ensg_uniprot_mapping: DataFrame mapping ENSG IDS to UniProt IDs.
    :type ensg_uniprot_mapping: pd.DataFrame
    :return: dictionary, with keys being tissues and values being a DataFrame of genes expressed in that tissue,
     annotated with GO-CC, and expression, specificity and secretome data
    :rtype: pd.DataFrame
    """
    logger.info('Annotating genes with expression and localisation features.')

    # gene_id (ENSG), uniprot. ensg_uniprot_mapping required because GO-CC uses UniProt, not ENSG
    ensembl_genes_xref = (ensembl_genes[['gene_id']]
                          .merge(ensg_uniprot_mapping.rename(columns={'ensg_id': 'gene_id'}), how='left'))

    # load HPA data onto genes
    logger.info('Loading HPA/GTEx expression and specificty data onto gene nodes.')
    ensembl_genes_features = (ensembl_genes_xref
                              .merge(tissue_exp
                                     .rename(columns={'Gene': 'gene_id', 'Source tissue': 'tissue'})
                                     .drop(columns='Gene name'), how='left')
                              .merge(tissue_specificity[['Gene', 'RNA tissue specificity', 'Source tissue']]
    .rename(
        columns={'Gene': 'gene_id', 'RNA tissue specificity': 'specificity_category', 'Source tissue': 'tissue'}),
        how='left')
                              .merge(secretome_data.rename(columns={'Gene': 'gene_id'}), how='left'))
    logger.info('Loading HPA/GTEx data completed successfully.')

    # add GO annotations
    logger.info('Loading GO-CC annotations onto gene nodes.')
    ensembl_genes_features = ensembl_genes_features.merge(
        go_cc_annots.rename(columns={'id': 'uniprot', 'go_id': 'go_cc'}), how='left').drop(columns='uniprot')
    logger.info('Loading GO-CC annotations completed successfully.')

    # fill in NAs
    ensembl_genes_features.fillna({
        'tissue': 'na',
        'specificity_category': '',
        'secretome_location': '',
        'secreted_function': '',
        'go_cc': ''
    }, inplace=True)

    # if one ENSG ID maps to many UniProts, there may be duplicate rows, each with the same ENSG ID but different GO IDs
    # collapse these duplicates down into one row and concat the GO IDs: one ENSG -> one set of GO IDs
    ensembl_genes_features_collapsed = (
        ensembl_genes_features.groupby([col for col in ensembl_genes_features.columns if col != 'go_cc'],
                                       as_index=False)
        .agg({'go_cc': lambda x: '|'.join(sorted(set(v for v in x if v)))})
    ).drop_duplicates()

    # change the 'tissue' to the format in the file names of the GTEx eQTLs
    logger.info('Harmonising how tissues in GTEx data are formatted.')
    ensembl_genes_features_collapsed['tissue'] = [re.sub(r'[^A-Za-z0-9]+', '_', x).strip('_')
                                                  if x != 'Brain - Spinal cord (cervical c-1)'
                                                  else 'Brain_Spinal_cord_cervical_c-1'
                                                  for x in ensembl_genes_features_collapsed['tissue'].to_list()]
    logger.info('Harmonisation complete.')

    # segregate by tissue
    ensembl_genes_by_tissue = {tissue: df for tissue, df in ensembl_genes_features_collapsed.groupby('tissue')}

    logger.info(
        'All data loaded onto gene nodes successfully, and gene nodes segregated by tissue-specific expression data.')

    # returns dict {tissue: df with gene_id (ENSG), go_cc, nTPM, tissue, specificity_category, secretome_location, secretome_function}
    return ensembl_genes_by_tissue


def stringdb_ensp_to_ensg(stringdb_df: pd.DataFrame) -> pd.DataFrame:
    """
    Map ENSPs found in StringDB to ENSGs.
    :param stringdb_df: DataFrame of StringDB interactions, with ENSP IDs
    :type stringdb_df: pd.DataFrame
    :return: Dataframe of StringDB interactions, with ENSGs instead of the original ENSPs
    :rtype: pd.DataFrame
    """
    # ENSP-ENSG mapping from STRING
    ensp_to_ensg = pd.read_csv('~/data/external/stringdb/9606.protein.aliases.v12.0.txt', sep='\t')
    ensp_to_ensg.columns = ['string_id', 'alias', 'source']
    ensp_to_ensg = ensp_to_ensg[ensp_to_ensg['alias'].str.contains('ENSG', na=False)][
        ['string_id', 'alias']].drop_duplicates(
        ignore_index=True)

    # convert STRING ENSP-based IDs to ENSGs
    stringdb_df_mapped = (stringdb_df
                          .merge(ensp_to_ensg
                                 .rename(columns={'string_id': 'protein1'}))
                          .drop(columns='protein1').rename(columns={'alias': 'gene1'})
                          .merge(ensp_to_ensg
                                 .rename(columns={'string_id': 'protein2'}))
                          .drop(columns='protein2').rename(columns={'alias': 'gene2'})).drop_duplicates()

    return stringdb_df_mapped


########################################################################################################################
#
# BUILD GRAPH
#
########################################################################################################################


def graph_stats_calculator(directed_graph: nx.MultiDiGraph) -> Dict[str, int]:
    """
    Calculate basic graph-wide statistics.
    :param directed_graph: Heterogeneous graph built with the data read in and processed in this script.
    :type directed_graph: nx.MultiDiGraph
    :return: Dictionary of graph statistics
    :rtype: Dict[str, int]
    """
    logger.info('Calculating basic stats on a nx.MultiDiGraph.')

    directed_graph_undirected = directed_graph.to_undirected()

    # homogenous gene graph for clustering - project
    gene_subgraph = directed_graph_undirected.subgraph(
        [n for n, d in directed_graph_undirected.nodes(data=True) if d.get("type") == 'gene'])
    gene_subgraph_simple = nx.Graph(gene_subgraph)

    largest_connected_component = directed_graph_undirected.subgraph(
        max(nx.connected_components(directed_graph_undirected), key=len))
    n_nodes = directed_graph.number_of_nodes()
    stats = {
        'n_nodes': n_nodes,
        'n_edges': directed_graph.number_of_edges(),
        'density': nx.density(directed_graph),
        'avg_degree': sum(dict(directed_graph.degree()).values()) / n_nodes,
        'max_degree': max(dict(directed_graph.degree()).values()),
        'n_components': nx.number_connected_components(directed_graph_undirected),
        'largest_component_size': max(len(x) for x in nx.connected_components(directed_graph_undirected)),
        'avg_clustering_coeff': nx.average_clustering(gene_subgraph_simple),
        'transitivity': nx.transitivity(gene_subgraph_simple)
    }

    logger.info('Stats calculation completed successfully.')

    return stats


def construct_single_tissue_graph(tissue: str, genes_by_tissue: pd.DataFrame, variant_data: pd.DataFrame,
                                  disease_nodes: pd.DataFrame, gwas_edges: pd.DataFrame, magma_edges: pd.DataFrame,
                                  proteomics_edges: Dict[str, pd.DataFrame], pqtl_edges: pd.DataFrame,
                                  eqtl_paths: List[str], interactions: pd.DataFrame,
                                  variant_gene_proximity_edges: pd.DataFrame, coexpressions: pd.DataFrame) \
        -> Dict[str, Union[nx.MultiDiGraph, Dict[str, int], pd.DataFrame]]:
    """
    Construct a networkx MultiDiGraph from tissue-specific data collated in this pipeline.
     Calculate basic graph stats on the constructed graph.
    :param tissue: Which tissue the graph is being constructed for
    :type tissue: str
    :param genes_by_tissue: Annotated Ensembl genes that are expressed in tissue
    :type genes_by_tissue: pd.DataFrame
    :param variant_data: variant hit data in hg38 coordinates
    :type variant_data: pd.DataFrame
    :param disease_nodes: DataFrame of ARD information
    :type disease_nodes: pd.DataFrame
    :param gwas_edges: DataFrame mapping Variant -> Disease GWAS links
    :type gwas_edges: pd.DataFrame
    :param magma_edges: DataFrame mapping Gene -> Disease MAGMA links
    :type magma_edges: pd.DataFrame
    :param proteomics_edges: DataFrame mapping Gene -> Disease proteomics links. [0]: Firth. [1]: Cox.
    :type proteomics_edges: Dict[str, pd.DataFrame]
    :param pqtl_edges: DataFrame of variants that are pQTLs, their pGenes and their effect sizes on the pGenes
    :type pqtl_edges: pd.DataFrame
    :param eqtl_paths: list of paths pointing to files containing tissue-specific GTEx eQTL data
    :type eqtl_paths: List[str]
    :param interactions: DataFrame of StringDB high-confidence (>=700) physical interactions, with ENSG IDs
    :type interactions: pd.DataFrame
    :param variant_gene_proximity_edges: DataFrame mapping variants to their closest genes, with distance
    :type variant_gene_proximity_edges: pd.DataFrame
    :param coexpressions: DataFrame of StringDB high-confidence (>=700) coexpression links, with ENSG IDs
    :type coexpressions: pd.DataFrame
    :return: Dictionary contaning (1) the constructed graph (2) basic statistics on the graph
     (3) DataFrame of graph nodes (4) DataFrame of graph edges
    :rtype: Dict[str, Union[nx.MultiDiGraph, Dict[str, int], pd.DataFrame]]
    """
    # everything is pinned on which proteins are expressed
    # eQTLs only from that tissue
    # everything else is tissue-invariant but will be subset depending on the proteins/eQTLs for that tissue

    logger.info(f"Constructing cellular network for {tissue}.")
    logger.info(f"Making data specific for {tissue}.")

    # eQTLs only from that tissue - unique_id, gene_id, effect
    eqtl_path = [x for x in eqtl_paths if tissue in x][0]
    eqtl_edges = map_variants_to_eqtls(variant_data=variant_data, eqtl_path=eqtl_path, eqtl_tissue=tissue)
    eqtl_genes = list(eqtl_edges['gene_id'].unique())

    # pQTLs
    pqtl_genes = list(pqtl_edges['ensg_id'].unique())

    # identify genes (1) expressed in this tissue (nTPM > 0.1) OR (2) where gene has nTPM < 0.1 but is e/pGene where
    #  e/pQTL has positive beta
    non_expressed_egenes_with_pos_beta = list(set(genes_by_tissue[genes_by_tissue['nTPM'] <= 0.1]['gene_id']) &
                                              set(eqtl_edges[eqtl_edges['effect'] > 0]['gene_id']))
    non_expressed_pgenes_with_pos_beta = list(set(genes_by_tissue[genes_by_tissue['nTPM'] <= 0.1]['gene_id']) &
                                              set(pqtl_edges[pqtl_edges['effect'] > 0]['ensg_id']))
    tissue_specific_genes = genes_by_tissue[(genes_by_tissue['nTPM'] > 0.1) |
                                            (genes_by_tissue['gene_id'].isin(non_expressed_egenes_with_pos_beta)) |
                                            (genes_by_tissue['gene_id'].isin(
                                                non_expressed_pgenes_with_pos_beta))].copy()
    tissue_specific_gene_list = list(tissue_specific_genes['gene_id'].to_list())

    # interactions, coexpressions
    tissue_specific_interactions = interactions[(interactions['gene1'].isin(tissue_specific_gene_list)) &
                                                (interactions['gene2'].isin(tissue_specific_gene_list))].copy()
    tissue_specific_coexpressions = coexpressions[(coexpressions['gene1'].isin(tissue_specific_gene_list)) &
                                                  (coexpressions['gene2'].isin(tissue_specific_gene_list))].copy()

    # variants proximal edges - unique_id, ensg_id
    tissue_specific_proximal_variants = variant_gene_proximity_edges[
        variant_gene_proximity_edges['gene_id'].isin(tissue_specific_gene_list)].copy()

    # variants - eQTLs (unique_id from variant_data), proximal variants (unique_id, ensg_id), pqtls (unique_id, etc.)
    tissue_specific_variant_ids = list(set(eqtl_edges['unique_id']) |
                                       set(tissue_specific_proximal_variants['unique_id']) |
                                       set(pqtl_edges['unique_id']))
    tissue_specific_variants = variant_data[
        variant_data['unique_id'].isin(tissue_specific_variant_ids)][['unique_id', 'MarkerID']].drop_duplicates().copy()

    # proteomics edges - ensg_id, disease, effect
    tissue_specific_proteomics_firth_edges = proteomics_edges['firth'][
        proteomics_edges['firth']['ensg_id'].isin(tissue_specific_gene_list)].copy()
    tissue_specific_proteomics_cox_edges = proteomics_edges['cox'][
        proteomics_edges['cox']['ensg_id'].isin(tissue_specific_gene_list)].copy()

    # MAGMA edges - ensg_id, disease
    tissue_specific_magma_edges = magma_edges[magma_edges['ensg_id'].isin(tissue_specific_gene_list)].copy()

    # GWAS edges - unique_id, disease, effect
    tissue_specific_gwas_edges = gwas_edges[gwas_edges['unique_id'].isin(tissue_specific_variant_ids)].copy()

    # pQTL edges
    tissue_specific_pqtl_edges = pqtl_edges[pqtl_edges['ensg_id'].isin(tissue_specific_gene_list)].copy()

    # eQTL edges
    tissue_specific_eqtl_edges = eqtl_edges[eqtl_edges['gene_id'].isin(tissue_specific_gene_list)].copy()

    # diseases - disease (ICD-10 three-letter code)
    tissue_specific_diseases = disease_nodes[
        disease_nodes['disease'].isin(
            list(set(tissue_specific_proteomics_firth_edges['disease']) |
                 set(tissue_specific_proteomics_cox_edges['disease']) |
                 set(tissue_specific_magma_edges['disease']) |
                 set(tissue_specific_gwas_edges['disease'])))].copy()

    logger.info(f"Making data specific for {tissue} completed successfully.")

    """ networkx object """
    # collate nodes
    logger.info('Collating node information.')
    nodes = pd.concat([tissue_specific_diseases.rename(columns={'disease': 'node_id'}).assign(node_type='disease'),
                       tissue_specific_variants.rename(columns={'unique_id': 'node_id'}).assign(node_type='variant'),
                       tissue_specific_genes.rename(columns={'gene_id': 'node_id'}).assign(node_type='gene')])
    # prefix features with 'features_'
    for node_column in nodes.columns:
        if node_column not in ['node_id', 'node_type']:
            nodes.rename(columns={node_column: 'feature_' + node_column}, inplace=True)
    # rename node_type to type
    nodes.rename(columns={'node_type': 'type'}, inplace=True)
    # check there are no duplicate nodes
    n_duplicated_nodes = nodes.duplicated(keep=False).sum()
    if n_duplicated_nodes > 0:
        logger.error(f"{n_duplicated_nodes} duplicated nodes found.")
        raise ValueError(f"{n_duplicated_nodes} duplicated nodes found.")
    # generate integer networkx IDs for each node
    # nodes['nx_id'] = np.arange(len(nodes))
    # generate maps between networkx IDs and human-readable IDs
    # nodes_map_nx_id_to_node_id = dict(zip(nodes['nx_id'], nodes['node_id']))
    # nodes_map_node_id_to_nx_id = dict(zip(nodes['node_id'], nodes['nx_id']))
    logger.info('Node information collated successfully.')

    # collate edges
    logger.info('Collating edge information.')
    edges = pd.concat([tissue_specific_eqtl_edges.rename(columns={'unique_id': 'source', 'gene_id': 'target'})
                      .assign(source_type='variant', target_type='gene', edge_type='eqtl_for'),
                       tissue_specific_pqtl_edges.rename(columns={'unique_id': 'source', 'ensg_id': 'target'})
                      .assign(source_type='variant', target_type='gene', edge_type='pqtl_for'),
                       pd.concat([tissue_specific_interactions,
                                  tissue_specific_interactions.assign(gene1=tissue_specific_interactions['gene2'],
                                                                      gene2=tissue_specific_interactions['gene1'])])
                      .rename(columns={'gene1': 'source', 'gene2': 'target'})
                      .assign(source_type='gene', target_type='gene', edge_type='interacts_with'),
                       pd.concat([tissue_specific_coexpressions,
                                  tissue_specific_coexpressions.assign(gene1=tissue_specific_coexpressions['gene2'],
                                                                       gene2=tissue_specific_coexpressions['gene1'])])
                      .rename(columns={'gene1': 'source', 'gene2': 'target'})
                      .assign(source_type='gene', target_type='gene', edge_type='coexpresses_with'),
                       tissue_specific_proximal_variants.rename(columns={'unique_id': 'source', 'gene_id': 'target'})
                      .assign(source_type='variant', target_type='gene', edge_type='proximal_to'),
                       tissue_specific_proteomics_firth_edges.rename(columns={'ensg_id': 'source', 'disease': 'target'})
                      .assign(source_type='gene', target_type='disease', edge_type='proteomics_firth_associated_with'),
                       tissue_specific_proteomics_cox_edges.rename(columns={'ensg_id': 'source', 'disease': 'target'})
                      .assign(source_type='gene', target_type='disease', edge_type='proteomics_cox_associated_with'),
                       tissue_specific_magma_edges.rename(columns={'ensg_id': 'source', 'disease': 'target'})
                      .assign(source_type='gene', target_type='disease', edge_type='magma_associated_with'),
                       tissue_specific_gwas_edges.rename(columns={'unique_id': 'source', 'disease': 'target'})
                      .assign(source_type='variant', target_type='disease', edge_type='gwas_associated_with')])
    # prefix features with 'features_'
    for edge_column in edges.columns:
        if edge_column not in ['source', 'target', 'source_type', 'target_type', 'edge_type']:
            edges.rename(columns={edge_column: 'feature_' + edge_column}, inplace=True)

    # check there are no duplicate edges
    duplicated_edges = edges[edges.duplicated(subset=['source', 'target', 'edge_type'])]
    duplicated_edges_unintentional = duplicated_edges[
        ~duplicated_edges['edge_type'].isin(['interacts_with', 'coexpresses_with'])]
    n_duplicated_edges = len(duplicated_edges_unintentional)
    if n_duplicated_edges > 0:
        logger.error(f"{n_duplicated_edges} duplicated edges found.")
        nodes.to_csv(HOME_APOCRITA_DIR + 'cellular_network/nodes_debugging.csv', index=False)
        edges.to_csv(HOME_APOCRITA_DIR + 'cellular_network/edges_debugging.csv', index=False)
        raise ValueError(f"{n_duplicated_edges} duplicated edges found.")
    # success
    logger.info('Edge information collated successfully.')

    # initialise networkx graph from edge dataframe
    logger.info('Initialising nx.MultiDiGraph.')
    single_tissue_graph = nx.from_pandas_edgelist(
        df=edges,
        source='source',
        target='target',
        edge_attr=[x for x in edges.columns if x not in ['source', 'target']],
        create_using=nx.MultiDiGraph()
    )

    # map nx_id to node attributes
    node_mapping = nodes.set_index('node_id')[[x for x in nodes.columns if x != 'node_id']].to_dict(orient="index")
    # add unconnected nodes
    single_tissue_graph.add_nodes_from(node_mapping.keys())
    # add node attributes
    nx.set_node_attributes(single_tissue_graph, node_mapping)
    logger.info(f"Single-tissue cellular network for {tissue} constructed successfully.")

    # basic stats on the graph
    stats = graph_stats_calculator(single_tissue_graph)

    # return built graph
    return {
        'graph': single_tissue_graph,
        'stats': stats,
        'nodes': nodes,
        'edges': edges
    }


########################################################################################################################
#
# MAIN EXECUTION
#
########################################################################################################################


def build_function():
    """
    Execute the graph build and export results.
    :return:
    """
    """ Function dumps """
    # read in data
    data_dict = read_in_data()
    # pull current Ensembl gene set
    # THESIS USED v115
    ensembl_gene_set = pull_current_ensembl_gene_set(current_ensembl_version=data_dict['current_ensembl_version'])
    # generate general UniProt <-> ENSG mappings, for all ENSGs
    uniprot_ensg_mapping = generate_ensg_uniprot_mappings(current_ensembl_version=data_dict['current_ensembl_version'],
                                                          ensembl_gene_set=ensembl_gene_set)
    # construct UniProt -> ENSG mappings for the Olink proteins: uniprot, ensg_id, gene_symbol, protein_field
    olink_uniprot_ensg_mappings = id_mappings_for_olink_proteins(uniprot_xref=uniprot_ensg_mapping)
    # map the GWAS hits from hg19 to hg38: imported hit data with pos column now hg38 and a unique_id column
    hits_with_hg38 = map_hg19_to_hg38(data_dict['hits'], chr_col='CHR', pos_col='POS')
    # link GWAS hits to disease: unique_id, effect, disease
    gwas_edges = make_gwas_edges(hits_with_hg38)
    # link genes to disease via MAGMA analyses
    magma_edges = make_magma_edges(current_ensembl_version=data_dict['current_ensembl_version'],
                                   magma_df=data_dict['magma'])
    # link proteins to disease via Firth and Cox analyses: ensg_id, disease, effect
    firth_edges = make_proteomics_edges(data_dict['firth'], olink_uniprot_ensg_mappings, algo='firth')
    cox_edges = make_proteomics_edges(data_dict['cox'], olink_uniprot_ensg_mappings, algo='cox')
    # link GWAS hits to genes via proximity: unique_id, gene_id
    variant_in_gene_edges = variant_gene_proximity_mapping(hits_with_hg38, ensembl_gene_set)
    # identify GWAS hits that are pQTLs: unique_id, ensg_id, effect
    pqtls = map_variants_to_pqtls(variant_data=hits_with_hg38, pqtl_data=data_dict['pqtls'],
                                  olink_uniprot_ensg_mapping=olink_uniprot_ensg_mappings)
    # create disease nodes (just the ICD10 three-letter code)
    disease_nodes = data_dict['disease_info'][['icd10_three_letter']].rename(columns={'icd10_three_letter': 'disease'})
    # create annotated gene nodes grouped in a dict by GTEx tissue: gene_id, go_cc, nTPM, tissue (expression),
    #  specificity_category (expression), secretome_location, secretome_function
    gene_nodes = annotate_genes(ensembl_genes=ensembl_gene_set,
                                go_cc_annots=data_dict['go_cc'],
                                tissue_exp=data_dict['hpa_exp'],
                                tissue_specificity=data_dict['hpa_tissue_specificity'],
                                secretome_data=data_dict['hpa_secretome'],
                                ensg_uniprot_mapping=uniprot_ensg_mapping)

    logger.info('Mapping ENSP IDs to ENSG IDs in StringDB annotations.')
    # pull high-confidence coexpression links from StringDB: gene1, gene2, coexpression
    coexpressions = stringdb_ensp_to_ensg(data_dict['stringdb_coexp'])
    # pull high-confidence interaction links from StringDB: gene1, gene2, experiments
    interactions = stringdb_ensp_to_ensg(data_dict['stringdb_int'])
    logger.info('Mapping ENSP IDs to ENSG IDs in StringDB completed successfully.')

    """ Build graphs for each tissue type """
    for tissue in data_dict['gtex_tissues']:
        if gene_nodes.get(tissue, None) is None:
            logger.info(f"There is no gene data for {tissue}.")
            continue
        graph_dict = construct_single_tissue_graph(tissue=tissue,
                                                   magma_edges=magma_edges,
                                                   proteomics_edges={'firth': firth_edges, 'cox': cox_edges},
                                                   genes_by_tissue=gene_nodes[tissue], variant_data=hits_with_hg38,
                                                   disease_nodes=disease_nodes, gwas_edges=gwas_edges, pqtl_edges=pqtls,
                                                   eqtl_paths=data_dict['eqtl_paths'], interactions=interactions,
                                                   variant_gene_proximity_edges=variant_in_gene_edges,
                                                   coexpressions=coexpressions)
        os.makedirs('~/data/internal/knowledge_graph/tissue_graphs/' + tissue + '/', exist_ok=True)
        nx.write_graphml(graph_dict['graph'], '~/data/internal/knowledge_graph/tissue_graphs/' + tissue + '/graph.graphml')
        (pd.DataFrame.from_dict(graph_dict['stats'], orient='index', columns=['value'])
         .assign(stat=graph_dict['stats'].keys())[['stat', 'value']]
         .to_csv('~/data/internal/knowledge_graph/tissue_graphs/' + tissue + '/graph_stats.csv', index=False))
        graph_dict['nodes'].to_csv('~/data/internal/knowledge_graph/tissue_graphs/' + tissue + '/nodes.csv', index=False)
        graph_dict['edges'].to_csv('~/data/internal/knowledge_graph/tissue_graphs/' + tissue + '/edges.csv', index=False)
        logger.info(f"Single-tissue cellular network for {tissue} and its stats, nodes and edges exported successfully "
                    f"to {'~/data/internal/knowledge_graph/tissue_graphs/' + tissue}.")


if __name__ == '__main__':
    build_function()
