"""
Test file for create_unfolded_sfs_from_df
"""


import pytest
from pandas import DataFrame
from .sfs_utils import create_unfolded_sfs_from_df


def test_create_unfolded_sfs_from_df_empty_df():
    input_df = DataFrame()  # Create an empty list

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="DataFrame is empty"):
        create_unfolded_sfs_from_df(input_df, "altcount", 10)


def test_create_unfolded_sfs_from_df_test_case():
    input_df = DataFrame({
        "chrom": "chr2L",
        "pos": [467, 891, 4427, 5192, 5391, 5579, 5735, 6266, 7271, 8323],
        "refcount": [5, 4, 1, 7, 5, 2, 6, 8, 9, 8],
        "altcount": [5, 6, 9, 3, 5, 8, 4, 2, 1, 2],
    })

    input_sum = 10.0

    expected_output = [0, 1, 2, 1, 1, 2, 1, 0, 1, 1, 0]
    expected_output_sum = sum(expected_output)

    # Assert that the unfolded SFS is a list
    assert create_unfolded_sfs_from_df(input_df, "altcount", 10) == expected_output
    assert input_sum == expected_output_sum
