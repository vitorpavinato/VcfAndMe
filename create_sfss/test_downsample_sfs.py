"""
Test file for downsample_sfs
"""


import pytest
from .sfs_utils import downsample_sfs

# Total number of gene copies in the original sfs
ORIGINAL_SIZE_SFS = 11

# Total number of gene compies in the downsampled sfs
DOWNSAMPLE_SIZE_SFS = 10


def test_downsample_sfs_empty_input_list():
    input_data = []  # Create an list

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="SFS is empty"):
        downsample_sfs(input_data, ORIGINAL_SIZE_SFS, DOWNSAMPLE_SIZE_SFS)


def test_downsample_sfs_test_sfs():
    input_sfs = [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    input_sum = 5.0
    expected_output_sfs = [0.09090909090909091, 1.090909090909091,
                           1.0909090909090908, 1.0909090909090908,
                           1.0909090909090908, 0.5454545454545454,
                           0.0, 0.0, 0.0, 0.0, 0.0]
    expected_output_sum = 5.0

    assert downsample_sfs(input_sfs,
                          ORIGINAL_SIZE_SFS,
                          DOWNSAMPLE_SIZE_SFS) == expected_output_sfs

    assert input_sum == expected_output_sum
