"""
Auxiliary functions for vcf_ann_to_table.py
"""

import os
import re
from pathlib import Path
from typing import List


# Check if the input file
def check_input_file(inputfile: str) -> str:
    """
    Simple function to check:
    - if the input file was provided by the user and if it exists;
    - if the input file has something other than the header; 
    - if the input file has at least 3 elements in the INFO fields;
    - the 3 elements in the INFO field should be::
        - AA, AC, and AF
    """

    # Check if input files are provided
    if inputfile == "" or inputfile is None:
        raise ValueError("file must be provided")

    # Check if the file exists
    if not os.path.exists(inputfile):
        raise ValueError("file does not exist")

    # Process the first input line
    first_line = None

    # Check if the input file has the right INFO fields format
    with open(inputfile, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if not line.startswith("##") and not line.startswith("#"):
                first_line = line
                break

    if first_line is not None:
        fields = first_line.strip().split('\t')
        # Check the INFO fields
        info_field = fields[7].split(';')

        # Check if the INFO fields:
        # Should have at least 3 elements in the list: AA, AC and AF
        if len(info_field) >= 3:
            if not (any(re.search(r'AA=.*', element) for element in info_field) and
                    any(re.search(r'AC=.*', element) for element in info_field) and
                    any(re.search(r'AF=.*', element) for element in info_field)):
                raise ValueError("Input is not supported by this script! It should have at least AA, AC and AF")

        else:
            # Handle the case when the INFO field has not enough elements
            raise ValueError("Input file has not enough information in the INFO field...")

    else:
        # Handle the case when the file is empty
        raise ValueError("Input file has only header information...")

    return inputfile


# Get the output file name
def output_file_name(inputfile: str, outputfile: str = None) -> str:
    """
    This function checks if the output file name is provided.
    If not, it will use the input file name with the .tsv extension
    """

    if outputfile == "" or outputfile is None:
        # Get the path and filename
        path, filename = os.path.split(inputfile)

        # Find the last occurrence of "/" in the path
        last_slash_index = path.rfind("/")

        # Extract the part before the last "/"
        base_path = path[:last_slash_index]

        # Append "tables" to the base path: the output file will be in the "tables" folder
        outputfile_path = base_path + "/tables"

        # Check if the outputfile path exists; create it if it doesn't
        if not Path(outputfile_path).exists():
            Path(outputfile_path).mkdir(parents=True)

        # Define the basename for the outputs from the filename
        basename, _ = filename.strip().split(".vcf")

        # Define the output file with a path
        outputfile = outputfile_path + "/" + basename + "_table.tsv"

    else:
        # Get the path and filename
        outputfile_path, filename = os.path.split(outputfile)

        # Check if the outputfile path exists; create it if it doesn't
        if not Path(outputfile_path).exists():
            Path(outputfile_path).mkdir(parents=True)

    return outputfile


# Write the processed VCF to a .TSV file
def write_tsv_file(lines: List[List[str]], header: List[str], outputfile: str) -> None:
    """
    Write the processed VCF to a .TSV file
    """

    # Prompt message
    print("Exporting the processed VCF to a .TSV file")
    with open(outputfile, "a", encoding="utf-8") as fo:
        fo.write(header + "\n")
        for line in lines:
            fo.write(("\t".join(str(item) for item in line)) + "\n")
    print("tsv file exported to: " + outputfile)
