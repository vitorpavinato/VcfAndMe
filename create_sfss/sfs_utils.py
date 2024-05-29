"""
This module contains a set of functions to create 
site-frequency spectrums (SFSs) and manipulate SFSs.
"""

import numpy as np
from pandas import DataFrame
from scipy.stats import hypergeom


# Create an unfolded SFS from a DataFrame
def create_unfolded_sfs_from_df(
    input_df: DataFrame,
    allele_count_column_name: str,
    max_number_of_haplotypes: int
) -> list[int | float]:
    """
    Create an unfolded SFS from a DataFrame.
    A DataFrame should have the at least the following columns:
    refcount and altcount.
    """

    # Raise error for an empty DataFrame
    if input_df.empty:
        raise ValueError("DataFrame is empty")

    # Create an empty unfolded SFS
    unfolded_sfs = [0] * (max_number_of_haplotypes + 1)

    # Fill up the unfolded SFS
    for _, row in input_df.iterrows():
        allele_counts = row[allele_count_column_name]
        unfolded_sfs[allele_counts] += 1

    return unfolded_sfs


# Fold an unfolded SFS
def fold_sfs(
    sfs: list[int | float]
) -> list[int | float]:
    """
    Fold an unfolded SFS
    """

    # Raise error for an empty SFS
    if not sfs:
        raise ValueError("SFS is empty")

    # Get the middle index
    mid_idx = int(len(sfs)/2)

    if len(sfs) % 2 > 0:
        # If the SFS is odd, get the middle value
        mid_value = sfs.pop(mid_idx)
    else:
        # If the SFS is even, set the middle value to None
        mid_value = None

    # Devide the SFS in two parts:
    # First part before the mid-index, the second part
    # after the mid-index. The two list should no include
    # the mid-index. Then revert the second list.
    before_mididx = sfs[0:mid_idx]
    after_mididx = sfs[mid_idx:]

    rev_after_mididx = after_mididx.copy()
    rev_after_mididx.reverse()

    fsfs = [x + y for x, y in zip(before_mididx, rev_after_mididx)]

    if mid_value is not None:
        fsfs.append(mid_value)

    return fsfs


# Get the a SFS from each of SNP total counts in dict
def create_sfs_dict(
    snp_counts_dict: dict[int, list[int]]
) -> dict[int, list[int]]:
    """
    Generate a SFS from each sample sizes (total counts) in dict.
    The SFS is a list of counts for each bin in the SFS. It expects
    a sorted dict of lists of altcounts for each total counts.
    It returns a sorted dict of SFS (which is a list of bin counts) for each
    total counts.
    """

    # Raise error for an empty dictionary
    if not snp_counts_dict:
        raise ValueError("Dictionary is empty")

    # Sort the dictionary of lists of altcounts for each total counts
    snp_counts_dict = dict(sorted(snp_counts_dict.items()))

    # Create an empty dict to save each popsize SFS
    # list of values for each key popsize
    snp_total_counts_dict = {}

    for key, value in snp_counts_dict.items():
        snp_total_counts_dict[key] = [0] * (key + 1)
        for alt_counts in value:
            snp_total_counts_dict[key][alt_counts] += 1

    return snp_total_counts_dict


# Downsampling SFSs
def downsample_sfs(
    original_sfs: list[int],
    original_size: int,
    sample_size: int
) -> list[int | float]:
    """
    Project a distribution an unfolded or folded site-frequency spectrum  
    to a sample size m < n.
    :param original_sfs: List of counts for each bin in the original SFS.
    Expects an SFS as a list of length sample size (n) + 1,
    with position i for i gene copies;
    :param popsize: Total number of gene copies in the original distribution.
    :param sampsize: Total number of gene copies in the
    downsampled distribution (m < n).
    :return: List of expected counts for each bin in the
    projected distribution.
    """

    # Raise error for an empty SFS
    if not original_sfs:
        raise ValueError("SFS is empty")

    # Create an empty sample SFS with sample_size + 1 bins
    sample_sfs = [0]*(sample_size+1)

    for pi, count in enumerate(original_sfs):
        for si in range(sample_size+1):
            prob = hypergeom.pmf(si, original_size, pi, sample_size)
            sample_sfs[si] += prob * count

    return sample_sfs


# Create a wrapper for SFS downsampling (and SFS folding) functions:
def downsample_sfs_in_dict(
    sfs_dict: dict,
    sample_size: int, fold: bool = False
) -> list[int | float]:
    """
    Project a distribution of an unfolded or folded site-frequency spectrum
    """

    # Raise error for an empty dictionary
    if not sfs_dict:
        raise ValueError("Dictionary is empty")

    # List of downsampled SFSs
    list_ds_sfs = []

    # Now donwsample only SFSs with key values higher than sample_size
    for key, values in sfs_dict.items():
        if key == sample_size:
            list_ds_sfs.append(values)
        else:
            ds_sfs = downsample_sfs(values, key, sample_size)
            list_ds_sfs.append(ds_sfs)

    # Conver the list of SFS to an np.array
    sfs_array = np.array(list_ds_sfs)

    # Sum column-wise to get the final SFS
    sfs = list(np.sum(sfs_array, 0))

    if fold:
        return fold_sfs(sfs)
    return sfs

# Extra functions to manipulate SFSs in dictionaries
# These functions were created to manipulate SFSs of
# DGRP2 datasets on 16/Feb/2024.

# Function to define common dict keys among a list of dictionaries
def find_common_dict_keys(list_of_chrm_dicts: list[dict]) -> set[str]:
    """
    Find the common keys among a list of dictionaries.
    The common keys are the SNP total counts common to
    all chromosome dictionaries.
    """

    common_keys = set(list_of_chrm_dicts[0].keys())
    for d in list_of_chrm_dicts[1:]:
        common_keys = common_keys.intersection(d.keys())
    return common_keys


# Function to sum the values of the common keys in all dictionaries.
def sum_sfs_in_dict_by_keys(
    list_of_chrm_dicts: list[dict],
    common_keys: set[str]
) -> dict:
    """
    Sum the values of the common keys in all dictionaries.
    """

    combined_chrm_sfs_by_key = {}

    for key in common_keys:
        combined_chrm_sfs_by_key[key] = [sum(x) for x in zip(*[d[key] for d in list_of_chrm_dicts])]

    return combined_chrm_sfs_by_key


# Function to fold the SFSs in dictionary
def fold_sfss_in_dict(dict_of_unfolded_sfss: dict) -> dict:
    """
    Fold the SFSs in dictionary of unfolded SFSs.
    """

    folded_sfs_dict = {}
    for key, values in dict_of_unfolded_sfss.items():
        folded_sfs = fold_sfs(values)
        folded_sfs_dict[key] = folded_sfs

    return folded_sfs_dict
