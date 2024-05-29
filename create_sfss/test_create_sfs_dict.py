"""
Test file for create_sfs_dict
"""


import numpy as np
import pytest
from .sfs_utils import create_sfs_dict


def test_create_sfs_dict_empty_input_dict():
    input_data = {}  # Create an empty dictionary

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="Dictionary is empty"):
        create_sfs_dict(input_data)


def test_create_sfs_dict_single_input():
    np.random.seed(42)  # Setting seed for reproducibility
    input_data = {11: [1, 2, 3, 4, 5],
                  12: [1, 2, 3, 4, 5]
                  }

    # Expected ouput
    expected_output = {11: [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], 
                       12: [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
                       }

    # Use the create_snp_total_counts_dict function
    assert create_sfs_dict(input_data) == expected_output
