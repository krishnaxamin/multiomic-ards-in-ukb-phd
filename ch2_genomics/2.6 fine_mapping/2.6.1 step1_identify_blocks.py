"""
For each disease, identify regions that contain associated variants. The variants in each region are less than 500kb from their neighbours.
Then, identify regions which are shared between disease pairs. Sharing = continguous or overlapping. Combine these regions to get regions of associated regions for each disease pair.
To lighten the overall computation for r-value computation, any regions that are contiguous or overlapping are combined, irrespective of disease pair or algorithm.
The r-values resulting can then be subset for the relevant ranges for each disease pair.
Also get disease-specific, algo-agnostic regions (combine regions for each disease from different algorithms, for colocalisation between differnt algorithms.
"""
from pandas import read_csv, merge, concat, DataFrame
from itertools import combinations

import numpy as np

# set maximum bp gap
max_gap = 1000000
buffer_size = 1000

# read in chr coordinates
chr_coordinates = read_csv('~/data/external/genomics/chr_coordinates_hg19.csv')

# read in SAIGE minimally adjusted data
saige_minimal = read_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_minimal.csv', dtype={'CHR': str})
saige_minimal = saige_minimal[(saige_minimal['assoc_label'] == 'assoc') & (saige_minimal['info_label'] == 'high_info')].copy().reset_index(drop=True)
saige_minimal = saige_minimal[['CHR', 'POS', 'disease', 'UKB_ref', 'UKB_alt']].copy()
saige_minimal['uniqueID'] = saige_minimal['CHR'].astype(str) + '_' + saige_minimal['POS'].astype(str) + '_' + saige_minimal['UKB_ref'] + '_' + saige_minimal['UKB_alt']
saige_minimal['CHR'] = saige_minimal['CHR'].astype(str)
saige_minimal.drop(['UKB_ref', 'UKB_alt'], axis=1, inplace=True)

# read in SAIGE lifestyles-adjusted data
saige_lifestyles = read_csv('~/data/internal/genomics/pan-ukbb-eur_assoc_saige_lifestyles.csv', dtype={'CHR': str})
saige_lifestyles = saige_lifestyles[(saige_lifestyles['assoc_label'] == 'assoc') & (saige_lifestyles['info_label'] == 'high_info')].copy().reset_index(drop=True)
saige_lifestyles = saige_lifestyles[['CHR', 'POS', 'disease', 'UKB_ref', 'UKB_alt']].copy()
saige_lifestyles['uniqueID'] = saige_lifestyles['CHR'].astype(str) + '_' + saige_lifestyles['POS'].astype(str) + '_' + saige_lifestyles['UKB_ref'] + '_' + saige_lifestyles['UKB_alt']
saige_lifestyles['CHR'] = saige_lifestyles['CHR'].astype(str)
saige_lifestyles.drop(['UKB_ref', 'UKB_alt'], axis=1, inplace=True)


def order_list_of_lists(list_of_lists, col_idx_to_sort):
    """
    Given a list of lists, return a list of the same low-level lists,
    with list components ordered by a particular element in the low-level list.
    :param list_of_lists:
    :param col_idx_to_sort:
    :return sorted_list:
    """
    df = DataFrame(list_of_lists, columns=list(range(len(list_of_lists[0]))))
    df.sort_values(by=col_idx_to_sort, inplace=True)
    sorted_list = df.values.tolist()
    return sorted_list


def disease_blocks(input_data, adjustment):
    """
    Make blocks for each disease.
    Within each block exist assoc variants which are within 1Mb of their neighbouring assoc variant.
    Blocks are padded with 1kb on each end.
    Blocks are referred to as 'relevant regions'.
    Blocks are used downstream, and are exported for use in fine mapping.
    :param input_data:
    :return relevant_regions:
    """

    relevant_regions = {}
    for disease in input_data.sort_values(by='disease').disease.unique():
        disease_assoc = input_data[(input_data.disease == disease)].copy().sort_values(by=['CHR', 'POS']).reset_index(drop=True)
        disease_regions = {}
        region_counter = 0
        for chromo in disease_assoc.CHR.unique():
            disease_assoc_chromo = disease_assoc[disease_assoc.CHR == chromo].copy().reset_index(drop=True)
            # for all variants on one chromosome to be in one region, the gaps between all variants must be less than 500kb
            min_bound = disease_assoc_chromo.POS.min() - buffer_size
            if min_bound < 0:
                min_bound = 0
            for i in range(len(disease_assoc_chromo) - 1):
                pos1 = disease_assoc_chromo.at[i, 'POS']
                pos2 = disease_assoc_chromo.at[i + 1, 'POS']
                # gap identified that is over 500kb
                if pos1 + max_gap < pos2:
                    max_bound = pos1 + buffer_size  # set max_bound for this region
                    disease_regions[region_counter] = [chromo, min_bound,
                                                       max_bound]  # save chr, min_ and max_bounds for this region
                    min_bound = pos2 - buffer_size  # redefine the new min_bound for the next region
                    region_counter += 1  # record that we are now looking at the next region

            # record the final max_bound. This happens irrespective of whether we found a 1Mb gap amongst variants on this chromosome
            max_bound = disease_assoc_chromo.POS.max() + buffer_size  # set max_bound as the max SNP co-ord + 1kb
            chr_max_bound = list(chr_coordinates[chr_coordinates['chr'] == 'chr' + chromo]['length'])[
                0]  # get the final co-ord of the whole chromosome
            if max_bound > chr_max_bound:  # check if that max_bound is outside the chromosome range
                max_bound = chr_max_bound  # if it is, set max_bound as the final co-ord of the chromosome
            disease_regions[region_counter] = [chromo, min_bound,
                                               max_bound]  # save chr, min_ and max_bounds for this region
            region_counter += 1  # record that we are now looking at the next region
        relevant_regions[disease] = disease_regions  # save the dict of regions into the overall dict, keyed by disease

    # export relevant regions for each disease
    all_regions_list = []
    for disease in relevant_regions.keys():
        region_df = DataFrame().from_dict(relevant_regions[disease], orient='index')
        region_df.columns = ['chr', 'start', 'end']
        region_df['disease'] = disease
        all_regions_list.append(region_df)
    all_regions_df = concat(all_regions_list)

    all_regions_df.to_csv('~/ch2_genomics/2.6 fine_mapping/disease_ld_indep_blocks_' + adjustment + '.csv', index=False)

    return relevant_regions


def pairwise_blocks(input_data, adjustment, relevant_regions):
    """
    Combine overlapping blocks for each pair of diseases.
    Blocks are referred to as 'relevant regions'.
    The resultant blocks are exported for colocalisation analyses of pairs of diseases that have overlapping blocks.
    :param input_data:
    :param relevant_regions:
    :return all_combos_regions:
    """

    disease_combinations = list(combinations(list(input_data.sort_values(by='disease').disease.unique()), r=2))

    # make relevant regions for each pair of disease.
    # overlapping regions from each disease indicate areas of common variants. these regions are merged into one.
    all_combos_regions = {}  # dict to store regions where variants common to both diseases in a pair are found
    for combo in disease_combinations:
        disease1 = combo[0]
        disease2 = combo[1]

        # make relevant regions for disease1
        # disease1_assoc = assoc[(assoc.disease == disease1)].copy().sort_values(by=['CHR', 'POS']).reset_index(drop=True)
        disease1_regions = relevant_regions[disease1]

        # make relevant regions for disease2
        # disease2_assoc = assoc[(assoc.disease == disease2)].copy().sort_values(by=['CHR', 'POS']).reset_index(drop=True)
        disease2_regions = relevant_regions[disease2]

        # identify chromosomes that have variants associated with both disease1 and disease2
        common_chromosomes = np.intersect1d([x[0] for x in disease1_regions.values()],
                                            [x[0] for x in disease2_regions.values()]).tolist()

        combo_regions = {}
        combo_region_counter = 0
        # iterate through regions from disease1 that are on a common_chromosome
        for region1 in [x for x in disease1_regions.values() if x[0] in common_chromosomes]:
            # look for regions from disease2 that are on the same chromosome as the region from disease1 and iterate through, looking for regions overlapping with region1
            for region2 in [x for x in disease2_regions.values() if x[0] == region1[0]]:
                if (region1[1] <= region2[2]) & (region1[1] >= region2[
                    1]):  # if the start of region1 is less than the end of region2, then overlap, with region1 > region2
                    combo_region_min_bound = min(region1[1], region2[1])  # min bound of combined region
                    combo_region_max_bound = max(region1[2], region2[2])  # max bound of combined region
                    combo_regions[combo_region_counter] = [region1[0], combo_region_min_bound,
                                                           combo_region_max_bound]  # save region in same format as before
                    combo_region_counter += 1  # record that a region has been saved and we're moving onto the next one
                elif (region2[1] <= region1[2]) & (region2[1] >= region1[
                    1]):  # if the start of region2 is less than the end of region1, then overlap, with region2 > region1
                    combo_region_min_bound = min(region1[1], region2[1])
                    combo_region_max_bound = max(region1[2], region2[2])
                    combo_regions[combo_region_counter] = [region1[0], combo_region_min_bound, combo_region_max_bound]
                    combo_region_counter += 1
        if len(combo_regions) > 0:
            all_combos_regions['-'.join(combo)] = combo_regions

    # export these regions for each disease combo
    all_combos_list = []
    for combo in all_combos_regions.keys():
        combo_df = DataFrame().from_dict(all_combos_regions[combo], orient='index')
        combo_df.columns = ['chr', 'start', 'end']
        combo_df['combo'] = combo
        all_combos_list.append(combo_df)
    all_combos_df = concat(all_combos_list)

    all_combos_df.to_csv('~/ch2_genomics/2.13 colocalisations/combo_regions_' + adjustment + '.csv', index=False)

    return all_combos_regions


def all_blocks_one_algo(relevant_regions):
    """
    For a given algorithm, combine all blocks that overlap regardless of disease.
    e.g. a block for E11 may overlap with a block for N18 - these should be combined.
    Blocks are referred to as 'relevant regions'.
    :param relevant_regions:
    :return algo_regions:
    """

    # extract all regions from their dictionaries
    list_of_disease_regions = []
    for x in relevant_regions.values():
        for xx in x.values():
            list_of_disease_regions.append(xx)
    list_of_disease_regions.sort()
    # get sorted list of all chromosomes which have these regions
    list_of_disease_chromos = sorted(list(set([x[0] for x in list_of_disease_regions])))

    algo_regions = {}
    algo_regions_counter = 0
    for chromo in list_of_disease_chromos:
        chromo_regions = [x for x in list_of_disease_regions if x[0] == chromo]
        algo_regions[algo_regions_counter] = [chromo, chromo_regions[0][1]]
        end_list = []  # collects ends of regions within the current 'final region'. max(end_list) used to determine the end of the 'final region'
        for i in range(len(chromo_regions)):
            end = chromo_regions[i][2]
            end_list.append(end)
            if i == len(chromo_regions) - 1:  # if looking at the last region in the chromosome-specific list
                # append the final end co-ordinate to the waiting 'final region' list in the algo_regions dict
                # algo_regions[algo_regions_counter].append(chromo_regions[len(chromo_regions) - 1][2])
                algo_regions[algo_regions_counter].append(max(end_list))
                algo_regions_counter += 1  # record that we are moving on to the next region
            else:  # not looking at the last region in the chromosome-specific list
                if chromo_regions[i + 1][1] <= max(
                        end_list):  # if the next region overlaps with the current 'final region', go directly to the next iteration
                    continue
                algo_regions[algo_regions_counter].append(
                    max(end_list))  # record the end co-ord of the current 'final region' list in the algo_regions dict
                algo_regions_counter += 1  # record that we are moving on to the next region
                algo_regions[algo_regions_counter] = [chromo, chromo_regions[i + 1][
                    1]]  # intialise the next 'final region' list in the algo_regions dict
                end_list = []  # clean end_list, ready for the next 'final region'

    return algo_regions


def all_blocks(algo_regions_list):
    """
    Combine all overlapping blocks for each analysis type together,
     to create a set of blocks which encompass variants associated with every disease by every algorithm,
     where the blocks follow the 1Mb-separation pattern.
    e.g. data from minimally adjusted has given a region which overlaps with a region from lifestyle-adjusted data -
     these should be combined.
    Resultant blocks are exported for use in LD calculations - LD r values will be calculated for each region separately.
    :param algo_regions_list:
    :return all_blocks:
    """
    # get list of all regions [[chr, start, end'], ['chr, start, end], ...]
    list_of_disease_regions = []
    for regions_list in algo_regions_list:
        for x in regions_list.values():
            list_of_disease_regions.append(x)

    # get sorted list of all chromosomes which have these regions
    list_of_disease_chromos = sorted(list(set([x[0] for x in list_of_disease_regions])))

    all_blocks_dict = {}
    all_blocks_dict_counter = 0
    for chromo in list_of_disease_chromos:
        chromo_regions = order_list_of_lists([x for x in list_of_disease_regions if x[0] == chromo], 1)
        all_blocks_dict[all_blocks_dict_counter] = [chromo, chromo_regions[0][1]]
        end_list = []  # collects ends of regions within the current 'final region'. max(end_list) used to determine the end of the 'final region'
        for i in range(len(chromo_regions)):
            end = chromo_regions[i][2]
            end_list.append(end)
            if i == len(chromo_regions) - 1:  # if looking at the last region in the chromosome-specific list
                # append the final end co-ordinate to the waiting 'final region' list in the all_blocks_dict dict
                # all_blocks_dict[all_blocks_dict_counter].append(chromo_regions[len(chromo_regions) - 1][2])
                all_blocks_dict[all_blocks_dict_counter].append(max(end_list))
                all_blocks_dict_counter += 1  # record that we are moving on to the next region
            else:  # not looking at the last region in the chromosome-specific list
                if chromo_regions[i + 1][1] <= max(
                        end_list):  # if the next region overlaps with the current 'final region', go directly to the next iteration
                    continue
                all_blocks_dict[all_blocks_dict_counter].append(
                    max(end_list))  # record the end co-ord of the current 'final region' list in the all_blocks_dict dict
                all_blocks_dict_counter += 1  # record that we are moving on to the next region
                all_blocks_dict[all_blocks_dict_counter] = [chromo, chromo_regions[i + 1][
                    1]]  # intialise the next 'final region' list in the all_blocks_dict dict
                end_list = []  # clean end_list, ready for the next 'final region'

    # export relevant regions for each disease
    all_blocks_dict_df = DataFrame(list(all_blocks_dict.values()), columns=['chr', 'start', 'end'])
    all_blocks_dict_df.to_csv('~/ch2_genomics/2.6 fine_mapping/blocks_for_ld_calcs.csv', index=False)

    return all_blocks_dict


saige_minimal_disease_blocks = disease_blocks(saige_minimal, 'minimal')
saige_minimal_combo_disease_blocks = pairwise_blocks(saige_minimal, 'minimal', saige_minimal_disease_blocks)
saige_minimal_all_blocks = all_blocks_one_algo(saige_minimal_disease_blocks)  # doesn't export

saige_lifestyles_disease_blocks = disease_blocks(saige_lifestyles, 'lifestyles')
saige_lifestyles_combo_disease_blocks = pairwise_blocks(saige_lifestyles, 'lifestyles', saige_lifestyles_disease_blocks)
saige_lifestyles_all_blocks = all_blocks_one_algo(saige_lifestyles_disease_blocks)  # doesn't export

all_blocks = all_blocks([saige_minimal_all_blocks, saige_lifestyles_all_blocks])
