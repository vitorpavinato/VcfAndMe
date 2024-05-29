"""
Annotate SNPs with the ancestral allele by parsimony only
if the outgroup allele is one of the target species allele.
This script is part of a collection of scripts for working 
vcf files. I tried to implement a linux principle of having
a tool to do a well a specific job. This is for rooting SNPs
only.
"""

import sys
import argparse
import subprocess  # this is being used to run samtools
from pathlib import Path  # This is being used to check if the file exists
from typing import List
from remake_vcf import create_samtools_faix


def add_root_info_to_vcf_header(
    input_file: Path, output_file: Path
) -> None:
    """
    Edit a vcf header to add root information. 
    This function makes the header to a standard format 
    compatible with GATK.
    """

    # Process header lines
    header_lines = []
    with open(input_file, 'r', encoding='utf-8') as input_file_obj:
        # Process lines starting with "##" or "#"
        for line in input_file_obj:
            if line.startswith("##") or line.startswith("#"):
                header_lines.append(line)

    # Add new lines and modify existing lines
    modified_header_lines = header_lines.copy()  # Create a copy of the original header lines

    # Remove two lines but add then later in the file processing
    # That is the reason to save then as variable for later
    # at same time I am using pop() method.
    chrm_line_header = modified_header_lines.pop(5)
    ref_line_header = modified_header_lines.pop(5)
    colum_names_line = modified_header_lines.pop(5)

    # Add new line
    modified_header_lines.append("##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n")

    # Add lines removed above
    modified_header_lines.append(chrm_line_header)
    modified_header_lines.append(ref_line_header)
    modified_header_lines.append(colum_names_line)

    # Write the modified header lines to the output file or perform further manipulations
    # For example, writing to the output file:
    with open(output_file, 'w', encoding='utf-8') as output_file_obj:
        for line in modified_header_lines:
            output_file_obj.write(line)


def root_snps_by_parsimony(
    input_file: Path = None,
    output_file: Path = None,
    outgroup_file: Path = None,
    samtools_path: Path = None
) -> None:
    """
    Root SNPs in a vcf based on the outgroup sequence.
    It takes a vcf file and return a new vcf with the rooted SNPs.
    """

    # Check if all required arguments are provided
    # This will raise an error if any of the required arguments are None

    # This is for the input_file
    if input_file is None:
        raise ValueError("input_file must be provided")

    # This is for the output_file
    if output_file is None:
        raise ValueError("output_file must be provided")

    if outgroup_file is None:
        raise ValueError("outgroup_file must be provided")

    # This is for the samtools_path
    if samtools_path is None:
        raise ValueError("samtools_path must be provided")

    # Check the existence of a outgroup reference file faidx
    create_samtools_faix(outgroup_file, samtools_path)
    print("Faidx created for the outgroup reference...")

    # run add_root_info_to_vcf_header
    add_root_info_to_vcf_header(input_file, output_file)
    print("Header updated with AA info...")

    # Process remaining lines
    with open(input_file, "r", encoding="utf-8") as input_file_obj, open(output_file, "a", encoding="utf-8") as output_file_obj:
        for line in input_file_obj:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            # Skip empty lines or lines with less than 9 fields
            if len(fields) < 9:
                continue

            # Create a list of alleles
            alleles = []

            # First add the reference allele
            alleles.append(fields[3])

            # Then extend the list to include the alternative allele(s)
            alleles.extend(fields[4].split(","))

            # Call samtools faidx to check outgroup reference base
            outgroup_ref_base = subprocess.run([
                samtools_path,
                "faidx",
                outgroup_file,
                (fields[0] + ":" + fields[1] + "-" + fields[1]),
            ],
                                    capture_output=True,
                                    check=False)

            # Get the reference base as a variable
            _, outgroup_ref_base, *_ = str(outgroup_ref_base.stdout).strip().split("\\n")

            # This will make sure that masked genome will work
            outgroup_ref_base = outgroup_ref_base.upper()

            # How I dealt with polarization is based on Machado et al. 2020.
            # The script tag the rooting only if the outgroup has the
            # alignment OR or has the one of the target species allelels.
            # Otherwise, the target species reference allele is used.
            if outgroup_ref_base in alleles:
                aa_annotation = f"AA={outgroup_ref_base}"
                if outgroup_ref_base == fields[3]:
                    fields[2] = 'root_ref'
                else:
                    fields[2] = 'root_alt'
            else:
                aa_annotation = f"AA={fields[3]}"
                fields[2] = 'root_unknown'

            # Include the aa_annotation string int fields[7]
            fields[7] = aa_annotation + ';' + fields[7]
            
            # Write the new line
            reassembled_line = '\t'.join(fields)

            # Write the new line
            with open(output_file, 'a', encoding='utf-8') as output_file_obj:
                output_file_obj.write(reassembled_line + "\n")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser("python root_vcf_by_parsimony.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input vcf file", dest="input_file")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Output vcf file", dest="output_file")
    parser.add_argument("-r", "--outgroup_ref_file", type=str, required=True, help="Outgroup reference file", dest="outgroup_file")
    parser.add_argument("-s", "--samtools_path", type=str, required=True, help="Path to samtools", dest="samtools_path")
    return parser


def main(argv: List[str]) -> None:
    """
    Main function
    """
    # Parse command line arguments
    parser = parse_args()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define input and output files
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    outgroup_file = Path(args.outgroup_file)
    samtools_path = Path(args.samtools_path)

    # Run root_vcf_by_parsimony
    root_snps_by_parsimony(input_file, output_file, outgroup_file, samtools_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
