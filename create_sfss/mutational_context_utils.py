'''
A set of function to create context dependent site-frequency spectrums (SFSs).
SNPs in the SFS might be constrained by distance, by similar mutational context
(here defined as the target SNP mutation plus and minos 1 nucleoted: AA/TC for
example). The functions included here help identify the distance contraint,
by calculating all the distances among synonymous and non-synonymous SNPs for
example ("calculate_distances"). It allows also identify blocks of consecutive
SNPs that might createambiguous mutational context
("find_consecutive_positions"). The remaining functions are part of the
mutational context pairing of the original pipeline, kept here for the record.
'''

from typing import Dict, List, Union, Tuple, Callable, Collection
import inspect
import numpy as np
from pandas import DataFrame
from sfs_utils import downsample_sfs, fold_sfs


def calculate_distances(
    df: DataFrame, type1: str, type2: str
) -> Collection[int | float]:
    """
    Giving a dataFrame, and two mutation types,
    this function returns the absolute distance
    between the type1 and type2. The dataFrame
    should have "effect" and "pos" fields.
    """

    # Filter rows for the specified types
    type1_df = df[df["effect"] == type1]
    type2_df = df[df["effect"] == type2]

    # Initialize a list to store distances
    distances = []

    # Iterate through each entry of type1
    for _, row1 in type1_df.iterrows():
        pos1 = row1["pos"]

        # Iterate through each entry of type2
        for _, row2 in type2_df.iterrows():
            pos2 = row2["pos"]

            # Calculate the absolute difference in positions
            distance = abs(pos1 - pos2)

            # Append the distance to the list
            distances.append(distance)

    return distances


def find_consecutive_positions(positions: list[int]) -> list[bool]:
    """
    This function find blocks of consecutive SNPs and
    returns False if a SNPs is part of a block of SNPs
    and True otherwise. This is useful to find and remove
    positions that might make ambiguous mutational
    contexts. Because of the nature of each of mutational
    context were defined based on the reference and
    alternative alleles.
    """
    pos = positions

    # Initialize a list to store consecutive blocks
    consecutive_blocks = []

    # Initialize variables for the current block
    current_block = [pos[0]]

    # Iterate through the positions starting from the second element
    for i in range(1, len(pos)):
        # Check if the current position is consecutive to the previous one
        if pos[i] - pos[i - 1] == 1:
            current_block.append(pos[i])
        else:
            # If not consecutive, start a new block
            consecutive_blocks.append(current_block)
            current_block = [pos[i]]

    # Append the last block
    consecutive_blocks.append(current_block)

    positions_to_keep = []
    for block in consecutive_blocks:
        if len(block) == 1:
            positions_to_keep.append(True)
        else:
            for i in block:
                positions_to_keep.append(False)

    return positions_to_keep


def set_snp_sequence_category(df: DataFrame) -> list[str]:
    """
    This define the mutational context each SNP in a dataFrame.
    Create a column containing the sequence category for each SNP
    in a DataFrame. The df should have, "ref", "alt" and "altflank"
    fields.
    """
    list_refalleles = df['ref']
    list_altalleles = df['alt']
    list_altflanks = df['altflank']
    
    sequence_categories = [f"{altflank[3]}{ref}/{alt}{altflank[5]}" for ref, alt, altflank  in zip(list_refalleles, list_altalleles, list_altflanks)]
    return sequence_categories


def create_sequence_categories_dict() -> dict[int, str]:
    """
    This creates a look-up dictionary for the sequence
    category. The dictionary then needs to be provided to
    "create_snp_dict" or "create_snp_dict_wrapper" in order
    to make a SNP dataFrame a dictionary of SNP organized by
    their SNP category.
    """

    # Possible SNPs:
    snp_types = [["A", "C"],
                 ["A", "G"],
                 ["A", "T"],
                 ["C", "G"],
                 ["C", "T"],
                 ["G", "T"]]

    # List of DNA nucleotides
    nucleotides = ["A", "C", "G", "T"]

    # Create an empty dictionary to save
    # the SNPs categories
    sequence_categories_dict = {}

    # Start the first dictionary key
    key_number = 1

    for snp_pairs in snp_types:
        for j in nucleotides:
            for i in nucleotides:
                key = key_number
                value = [
                         i + snp_pairs[0] + "/" + snp_pairs[1] + j,
                         i + snp_pairs[1] + "/" + snp_pairs[0] + j
                         ]
                sequence_categories_dict[key] = value

                key_number += 1
                # # Set the reciprocal
                # key = key_number+1
                # value = i + snp_pairs[1] + "/" + snp_pairs[0] + j
                # sequence_categories_dict[key] = value
                # key_number += 2

    return sequence_categories_dict


def create_snp_dict(
    df: DataFrame,
    sequence_categories_dict: dict[int, str]
) -> dict[int, dict[int, list[int, int]]]:
    """
    Covert a given pandas DataFrame containing SNPs
    to a dictionary. It needs a sequence category dictionary
    created with "create_sequence_categories_dic". The df
    should have "refcount", "altcount" and "sequence_category"
    fields.
    """

    # Initialize the snps dictionary
    snps_dict = {}

    # Fill up the SNPs dictionary using the table data
    for _, row in df.iterrows():
        pos = row['pos']
        ref_counts = row['refcount']
        alt_counts = row['altcount']
        total_counts = ref_counts + alt_counts
        sequence_category = row['sequence_category']

        # Find the corresponding numeric key(s) from the lookup dictionary
        numeric_keys = [
            key for key, value in sequence_categories_dict.items()
            if sequence_category in value
            ]

        # Add the data to the selected_snps dictionary for each numeric key
        for numeric_key in numeric_keys:
            # numeric_key_int = int(numeric_key)
            if numeric_key not in snps_dict:
                snps_dict[numeric_key] = {}
            snps_dict[numeric_key][pos] = [alt_counts, total_counts]

    return snps_dict


def create_snp_dict_wrapper(
    df: DataFrame,
    sequence_categories_dict: dict[int, str],
    snp_class: str
) -> dict[int, dict[int, list[int, int]]]:
    """
    Wrapper for create_snp_dict. This wrapper allows to process different
    effect mutations: introns, non-synonymous and synonymous.
    It expect the "effect" field to be one of "INTRON",
    "NON_SYNONYMOUS_CODING", or "SYNONYMOUS_CODING".
    """

    if snp_class == "introns":
        introns_df = df[df['effect'] == "INTRON"]
        introns_dict = create_snp_dict(introns_df, sequence_categories_dict)
        return introns_dict

    elif snp_class == "exons":
        nonsyns_df = df[df['effect'] == "NON_SYNONYMOUS_CODING"]
        syns_df = df[df['effect'] == "SYNONYMOUS_CODING"]
        nonsyns_dict = create_snp_dict(nonsyns_df, sequence_categories_dict)
        syns_dict = create_snp_dict(syns_df, sequence_categories_dict)
        return nonsyns_dict, syns_dict

    else:
        raise ValueError("snp_class must be 'introns' or 'exons'")


def is_sorted(lst):
    """
    From Codeium: Check if list is sorted
    """
    return all(lst[i] <= lst[i+1] for i in range(len(lst)-1))


def find_a_snp(inpute_list: list[int]) -> tuple[int, int] | None:
    """
    This simply take the first item in sorted list as the closest
    item to the other list. The list of position is sorted,
    so the first item should be the closest item to the first
    other list.
    """
    # Check if the list is sorted
    if not is_sorted(inpute_list):
        sorted_list = sorted(inpute_list)
    else:
        sorted_list = inpute_list

    # Find the closest position
    closest_num_idx = 0
    closest_num = sorted_list[0]

    return closest_num, closest_num_idx


def find_a_closest_snp(inpute_list: list[int], k: int) -> tuple[int, int] | None:
    """
    Find the closest position to a target position k in
    a list of sorted positions.
    This functions is part of the find_close_snp_pairs.
    """

    # Check if the list is sorted
    if not is_sorted(inpute_list):
        sorted_list = sorted(inpute_list)
    else:
        sorted_list = inpute_list

    # Find the closest position
    closest_num_idx = 0
    closest_num = sorted_list[0]
    if len(sorted_list) > 1:
        for idx, num in enumerate(sorted_list):
            if abs(num - k) < abs(closest_num - k):
                closest_num = num
                closest_num_idx = idx
            if num > k:
                break

    return closest_num, closest_num_idx


def find_a_closest_snp_within_interval(
    inpute_list: list[int],
    k: int,
    interval: int
) -> tuple[int, int] | None:
    """
    Find the closest position to a target position k in a list of sorted 
    positions within a predefined interval.
    This functions is part of the find_close_snp_pairs.
    """

    # Check if the list is sorted
    if not is_sorted(inpute_list):
        sorted_list = sorted(inpute_list)
    else:
        sorted_list = inpute_list

    # Check if interval is provided
    if not interval:
        raise ValueError("interval must be provided")

    # Find the closest position
    closest_num_idx = 0
    closest_num = sorted_list[0]
    if len(sorted_list) > 1:
        for idx, num in enumerate(sorted_list):
            if abs(num - k) < abs(closest_num - k) and abs(num - k) <= interval:
                closest_num = num
                closest_num_idx = idx
            if num > k:
                break

    return closest_num, closest_num_idx


# NEED SOME WORK IN HERE
def find_closest_snp_pairs(
    dict1: Dict[int, Dict[int, List[Tuple[int, int]]]],
    dict2: Dict[int, Dict[int, List[Tuple[int, int]]]]
    # ,
    # find_closest_func: Union[Callable[[List[int]], Tuple[int, int]], Callable[[List[int], int], Tuple[int, int]], Callable[[List[int], int, int], Tuple[int, int]]]
) -> Tuple[Dict[int, Dict[int, List[Tuple[int, int]]]], Dict[int, Dict[int, List[Tuple[int, int]]]]]:
    """Find pairs of closest SNPs from two distinct
    dictionaries. The first dictionary must be the one
    with less SNPs (usually introns or any other neutral
    control)
    """

    # Check if both dictionaries are not empty
    if not dict1 or not dict2:
        raise ValueError("One or both dictionaries are empty!")

    # # Check if find_closest_func is provided
    # if find_closest_func is None:
    #     raise ValueError("find_closest_func must be provided!")

    # # Check if find_closest_func is callable
    # if not callable(find_closest_func):
    #     raise ValueError("find_closest_func must be callable!")

    # # Define find_snp function with k and interval arguments
    # def selected_func(func, imput_list=None, k=None, interval=None) -> Callable:
    #     def wrapper(*args, **kwargs):
    #         if 'imput_list' in kwargs:
    #             return func(*args, **kwargs, imput_list=imput_list)
    #         elif 'k' in kwargs and 'interval' in kwargs and 'imput_list' in kwargs:
    #             return func(*args, **kwargs)
    #         elif 'k' in kwargs:
    #             return func(*args, **kwargs, k=k)
    #         elif 'interval' in kwargs:
    #             return func(*args, **kwargs, interval=interval)
    #         else:
    #             return func(*args, **kwargs, k=k, interval=interval)
    #     return wrapper

    # # Define the interval based on the find_closest_func
    # if len(inspect.signature(find_closest_func).parameters) == 3:
    #     interval = int(find_closest_func.__annotations__.get("interval"))
    #     find_snp = selected_func(find_closest_func, interval)
    # else:
    #     find_snp = selected_func(find_closest_func)

    # Two empty dictionaries to save each SNP
    # in a pair found in both input dictionaries
    paired_snps_dict1 = {}
    paired_snps_dict2 = {}

    # Transverse each sequence category (seq_key) in the first
    # dictionary, finding the same seq_key in the second dictionary
    for seq_key, dict1_positions in dict1.items():
        if seq_key in dict2:
            dict2_positions = dict2[seq_key]
            list_dict2_positions = list(dict2_positions)

            # When there is a seq_key in both, transvers the list of
            # positions in the first, finding a closes SNPs in the second.
            # Remove the kept SNPs from the second dictionary list to
            # avoid repeting the same SNP. Leave the loop when the list of
            # SNPs in the second dictionary seq_key is exhausted.

            for pos in dict1_positions:
                # Find the closest position in the second dictionary
                # This part accepts a function as an argument
                closest_pos, closest_pos_idx = find_a_snp(list_dict2_positions)

                paired_snps_dict1.setdefault(seq_key, {})[pos] = dict1_positions[pos]
                paired_snps_dict2.setdefault(seq_key, {})[closest_pos] = dict2_positions[closest_pos]

                # Remove the selected position from the list
                # This works like a FIFO algorithm
                list_dict2_positions.pop(closest_pos_idx)
                if len(list_dict2_positions) == 0:
                    break

    return paired_snps_dict1, paired_snps_dict2


# Maybe remove this function
def create_unfolded_sfs_from_snp_dict(
    snp_dict: dict[int, dict[int, list[int, int]]],
    count_data_type: str,
    max_sample_size: int,
    min_sample_size: int = None,
    folded: bool = False,
) -> list[int | float]:
    """
    Get the SFS from each sample size in dict. The input
    dictionary should have the following hierarchy:
    sequence_category: position: [alt_count, total_count].
    """

    # Raise error for an empty dictionary
    if not snp_dict:
        raise ValueError("Dictionary is empty")

    # Check if max_sample_size are provided
    if not max_sample_size:
        raise ValueError("at least max_sample_size must be provided")

    # Check if and max_sample_size are integers
    if not isinstance(max_sample_size, int):
        raise ValueError("max_sample_size must be integers")

    # Railroad pattern to determine what to do based on the the count_data_type
    if count_data_type == "imputed":  # imputed data needs max sample size only
        print("Working with imputed data...")

        # Create an empty unfolded SFS
        unfolded_sfs = [0] * (max_sample_size + 1)

        # Fill in the unfolded SFS
        for seq_key in snp_dict:
            for pos in snp_dict[seq_key]:
                alt_count = snp_dict[seq_key][pos][0]
                unfolded_sfs[alt_count] += 1

        if folded:
            return fold_sfs(unfolded_sfs)
        return unfolded_sfs

    elif count_data_type == "downsampled":  # downsample needs a min and max sample sizes

        # Check if max_sample_size are provided
        if min_sample_size is None:
            raise ValueError("For downsampled data max_sample_size must be provided")

        # Check if and max_sample_size are integers
        if not isinstance(min_sample_size, int):
            raise ValueError("min_sample_size must also be integers")

        # Check if min_sample_size < max_sample_size
        if min_sample_size > max_sample_size:
            raise ValueError("min_sample_size must be less than or equal to max_sample_size")

        print("Working with downsampled data...")

        # Create an empty dict to save each popsize SFS
        # list of values for each key popsize
        sample_size_dict = {}

        for i in range(min_sample_size, max_sample_size + 1):
            sample_size_dict[i] = [0] * (i+1)

        # Sequence category -wise loop to fill the
        # dictionary of popsize SFS
        for seq_key in snp_dict:
            dict_pos = snp_dict[seq_key]
            for pos in dict_pos:
                snp_counts = dict_pos[pos]
                sample_size_dict[snp_counts[1]][snp_counts[0]] += 1

        # Convert the sample_size_dict to a list and apply
        # downsampling in the SFS for higher values than
        # min_sample_size
        unfolded_sfs_list = []
        for sample_size, usfs in sample_size_dict.items():
            if sample_size == min_sample_size:
                unfolded_sfs_list.append(usfs)
            else:
                ds_usfs = downsample_sfs(usfs, sample_size, min_sample_size)
                unfolded_sfs_list.append(ds_usfs)

        # Conver the list of SFS to an np.array
        unfolded_sfs_array = np.array(unfolded_sfs_list)

        # Sum column-wise to get the final SFS
        unfolded_sfs = list(np.sum(unfolded_sfs_array, 0))

        if folded:
            return fold_sfs(unfolded_sfs)
        return unfolded_sfs
    else:
        raise ValueError("count_data_type must be 'downsampled' or 'imputed'")
