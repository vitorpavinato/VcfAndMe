"""
vcf_to_tsv_func.py
This module contains functions and classes to use in the vcf_to_tsv script.
They are lower level functions, not for the user. They are used to check
file paths, set output files names, etc.
"""


import os
import subprocess
from pathlib import Path
from typing import List


# Check if the input file
def check_input_file(inputfile: str) -> str:
    """
    Simple function to check:
    - if the input file was provided by the user and if it exists;
    - if the input file has something in the INFO fields;
    - if the input file has at least 4 elements in the INFO field;
    - The INFO field should have: AA, AC, AF and EFF.
    """

    # Check if input files are provided
    if inputfile == "" or inputfile is None:
        raise ValueError("file must be provided")

    # Check if the file exists
    if not os.path.exists(inputfile):
        raise ValueError("file does not exist")

    # Check if the input file has the right INFO fields format
    with open(inputfile, "r", encoding="utf-8") as input_file:
        first_line = None
        for line in input_file:
            if not line.startswith("##") and not line.startswith("#"):
                first_line = line
                break

    if first_line is not None:
        fields = first_line.strip().split('\t')
        # Check the INFO fields
        info_field = fields[7].split(';')

        # Check if the INFO fields have at least 3 elements or if it is not empty:
        # INFO field should have: AC, AF and EFF.
        if len(info_field) < 4 or fields[7] == "":
            raise ValueError("Input file has missing elements in INFO field. Was the vcf annotated with SNPEff?")

        # Handle the case when the file is not empty
        # and has at least 3 elements in INFO field
        # Split the elements in info_field_list
        first = info_field[0].split('=')[0]
        second = info_field[1].split('=')[0]
        third = info_field[2].split('=')[0]
        fourth = info_field[3].split('=')[0]

        # DONT NEED THIS ANYMORE
        if (first != "AA") and (second != "AC") and (third != "AF") and (fourth != "EFF"):
            raise ValueError("Input is not supported by this script! It should have at least AA, AC, AF and EFF")

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
        outputfile = outputfile_path + "/" + basename + "_table.vcf"

    else:
        # Get the path and filename
        outputfile_path, filename = os.path.split(outputfile)

        # Check if the outputfile path exists; create it if it doesn't
        if not Path(outputfile_path).exists():
            Path(outputfile_path).mkdir(parents=True)

    return outputfile


# Check samtools path
def check_samtools_path(samtools_path: str) -> str:
    """
    Check if the samtools path is provided
    """

    # Check if samtools path is provided
    if samtools_path is None:
        raise ValueError("samtools_path must be provided")

    samtools = samtools_path

    return samtools


# Check if the reference file is provided
def check_reference_genome_file(reference: str) -> str:
    """
    Check if the reference file is provided
    """

    # Check if the reference file is provided
    if reference is None:
        raise ValueError("reference must be provided")

    return reference


# Check if faidx associated to a reference genome file in fasta exists
def create_faidx(reference: str, samtools: str) -> None:
    """
    Check if the reference file has a faidx.
    Otherwise, create it.
    """

    # Check if .fai associated to a reference genome file in fasta exists
    if not os.path.exists((reference + ".fai")):
        print("Creating .fai associated to a reference genome file in fasta")
        faidx = subprocess.run([samtools, "faidx", reference], capture_output=True, check=False)
    else:
        print(".fai associated to a reference genome file in fasta already exists")


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
