"""
Test file for fold_sfs
"""


import pytest
from .sfs_utils import fold_sfs


def test_fold_sfs_empty_sfs():
    input_data = []  # Create an empty list

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="SFS is empty"):
        fold_sfs(input_data)


def test_fold_sfs_even_number_sample_sfs():
    input_data = [14, 9, 2, 2, 1, 1, 1]  # Create SFS of size n+1, n=6
    input_sum = 30.0

    expected_output = [15, 10, 3, 2]
    expected_output_sum = sum(expected_output)
    # Assert that the expected output is the same as
    # calculated output for an SFS of size n+1
    assert fold_sfs(input_data) == expected_output
    assert input_sum == expected_output_sum


def test_fold_sfs_odd_number_sample_sfs():
    input_data = [12, 8, 3, 2, 2, 1, 1, 1]  # Create SFS of size n+1, n=7
    input_sum = 30.0

    expected_output = [13, 9, 4, 4]
    expected_output_sum = sum(expected_output)
    # Assert that the expected output is the same as
    # calculated output for an SFS of size n+1
    assert fold_sfs(input_data) == expected_output
    assert input_sum == expected_output_sum
