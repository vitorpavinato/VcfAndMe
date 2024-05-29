"""
Test file for create_snp_total_counts_dict
"""


import numpy as np
import pandas as pd
import pytest
from .snp_utils import create_snp_total_counts_dict


def test_create_snp_total_counts_dict_empty_dataframe():
    input_data = pd.DataFrame()  # Create an empty DataFrame

    # Assert that calling create_snp_total_counts_dict with an empty 
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="DataFrame is empty"):
        create_snp_total_counts_dict(input_data)


def test_create_snp_total_counts_dict_single_input():
    np.random.seed(42)  # Setting seed for reproducibility
    input_data = pd.DataFrame({
        "chrom": "chr2L",
        "pos": np.sort(np.random.randint(1, 10001, size=10)),
        "refcount": [(11-1), (11-2), (11-3), (11-4), (11-5),
                     (12-1), (12-2), (12-3), (12-4), (12-5)],
        "altcount": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "effect": np.random.choice(["intron", "non-synonymous"], size=10),
    })

    # Expected ouput
    expected_output = {
                        11: [1, 2, 3, 4, 5],
                        12: [1, 2, 3, 4, 5]
                        }

    # Sort the vcf file
    input_data = input_data.sort_values(by=["pos"])

    # Use the create_snp_total_counts_dict function
    assert create_snp_total_counts_dict(input_data) == expected_output
