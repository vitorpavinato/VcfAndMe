"""
Test file for filter_snps_by_interval
"""


import pytest
from pandas import DataFrame
from pandas.testing import assert_frame_equal
from .snp_utils import filter_snps_by_interval

empty_tsv_df = DataFrame()  # Create an empty DataFrame
empty_bed_df = DataFrame()  # Create an empty DataFrame

# Generate snp_df DataFrame
tsv_df = DataFrame({
    "chrom": "chr2L",
    "pos": [467, 891, 4427, 5192, 5391, 5579, 5735, 6266, 7271, 8323],
    "refcount": [5, 4, 1, 7, 5, 1, 4, 0, 9, 5],
    "altcount": [8, 0, 9, 2, 6, 3, 8, 2, 4, 2],
    "effect": ["intron", "intron", "intron", "intron", "intron",
               "non-synonymous", "intron", "non-synonymous",
               "non-synonymous", "non-synonymous"],
})

# Generate random data for bed DataFrame
bed_df = DataFrame({
    "chrom": "chr2L",
    "starts": [1000, 4000, 5000, 7000, 9000],
    "ends": [1500, 4300, 6999, 7100, 9500],
    "interval": ["intron1", "intron2", "intron3", "intron4", "intron5"],
})


def test_filter_snps_by_interval_empty_snp_df():
    # Assert that calling filter_snps_by_interval with an empty
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="tsv_df is empty"):
        filter_snps_by_interval(empty_tsv_df, bed_df)


def test_filter_snps_by_interval_empty_bed_df():
    # Assert that calling filter_snps_by_interval with an empty
    # DataFrame raises a ValueError
    with pytest.raises(ValueError, match="bed_df is empty"):
        filter_snps_by_interval(tsv_df, empty_bed_df)


def test_filter_snps_by_interval_input_df():
    expected_output_tsv_df = DataFrame({
        "chrom": "chr2L",
        "pos": [5192, 5391, 5579, 5735, 6266],
        "refcount": [7, 5, 1, 4, 0],
        "altcount": [2, 6, 3, 8, 2],
        "effect": ["intron", "intron", "non-synonymous",
                   "intron", "non-synonymous"],
    })

    filtered_tsv_df = filter_snps_by_interval(tsv_df, bed_df)
    filtered_tsv_df.reset_index(drop=True, inplace=True)
    assert_frame_equal(filtered_tsv_df, expected_output_tsv_df)
