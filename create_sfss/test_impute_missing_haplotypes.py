"""
Test file for impute_missing_haplotypes
"""


import pytest
from pandas import DataFrame
from pandas.testing import assert_frame_equal
from .snp_utils import impute_missing_haplotypes


def test_impute_missing_haplotypes_empty_df():
    input_df = DataFrame()  # Create an empty list

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="DataFrame is empty"):
        impute_missing_haplotypes(input_df, 10)


def test_impute_missing_haplotypes_test_df():
    input_df = DataFrame(
        {'refcount': [1, 6, 3, 4],
         'altcount': [8, 2, 7, 2],
         'totalcount': [9, 8, 10, 6]})

    expected_output = DataFrame(
        {'refcount': [1, 8, 3, 8],
         'altcount': [9, 2, 7, 2],
         'totalcount': [9, 8, 10, 6]})

    output_df = impute_missing_haplotypes(input_df, 10)
    assert_frame_equal(output_df, expected_output)
