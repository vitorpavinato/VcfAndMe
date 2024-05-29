"""
Test file for create_unfolded_sfs_from_snp_dict
"""

import numpy as np
import pytest
from .mutational_context_utils import create_unfolded_sfs_from_snp_dict


def test_create_unfolded_sfs_from_snp_dict_empty_input_dict():
    input_data = {}  # Create an empty dictionary

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="Dictionary is empty"):
        create_unfolded_sfs_from_snp_dict(input_data, "imputed", 10)

# Test case for create_unfolded_sfs_from_snp_dict with non empty input
input_data = {
    1: {1: [1, 10], 2: [2, 10], 3: [8, 10]},
    2: {4: [1, 10], 5: [2, 10], 6: [9, 10]},
    3: {7: [1, 10], 8: [2, 10], 9: [10, 10]},
}

def test_create_unfolded_sfs_from_snp_dict_unfolded_sfs():
    # Expected ouput
    expected_output = [0, 3, 3, 0, 0, 0, 0, 0, 1, 1, 1]

    # Use the create_snp_total_counts_dict function
    assert create_unfolded_sfs_from_snp_dict(input_data, "imputed", 10, folded=False) == expected_output


def test_create_unfolded_sfs_from_snp_dict_folded_sfs():
    # Expected ouput
    expected_output = [1, 4, 4, 0, 0, 0]

    # Use the create_snp_total_counts_dict function
    assert create_unfolded_sfs_from_snp_dict(input_data, "imputed", 10, folded=True) == expected_output
