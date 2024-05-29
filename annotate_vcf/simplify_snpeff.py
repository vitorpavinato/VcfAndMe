"""
Sript to simplify SNPEff annotation
For the moment it only takes the first term
It might be redundant as SNPEff has this option
But the idea is to have a more probabilistic approach
in selecting annotated terms.
"""

import os
import argparse
import sys


# Function to check if a file exists
def parser_and_checker(file: str, outfile: str = None) -> bool:
    """
    Simple function to check:
    - if the input file was provided by the user and if it exists;
    - if the input file has the right INFO fields format.
    - if the input file has at least 4 elements in the INFO field
    related variables.
    """

    # Check if input files are provided
    if file == "" or file is None:
        raise ValueError("file must be provided")

    # Check if the file exists
    if not os.path.exists(file):
        raise ValueError("file does not exist")

    # Get the path and filename
    path, filename = os.path.split(file)

    # Define the basename for the outputs from the filename
    basename, _ = filename.strip().split('.vcf')

    # Define input and output files
    inputfile = file

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

        # Check if the INFO fields have at least 4 elements or if it is not empty:
        # INFO field should have: AA, AC, AF and EFF.
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
            raise ValueError("Input is not supported by this script!")

    else:
        # Handle the case when the file is empty
        raise ValueError("Input file has only header information...")

    # This handles the output file name
    if outfile is None:
        outputfile = path + "/" + basename + "_simplified.vcf"
    else:
        outputfile = outfile

    return inputfile, outputfile


# simplify_snpeff_default()
def simplify_snpeff_default(file: str, outfile: str = None,
                            keeponlyterms: bool = False,
                            method: str = "First") -> None:

    """
    This function simplifies SNPEff entries of a VCF file.
    This is the default method without custom annotations.
    For the moment, method don't do anything.
    """

    # Define input and output files
    inputfile, outputfile = parser_and_checker(file, outfile)

    # This handles if keeponlyterms is set to True
    if keeponlyterms:
        # Define annotation effects terms we care about
        relevant_effect_terms = ['INTRON', 'SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING', 'INTERGENIC']

    # Copy the header from the inputfile to the ouputfile
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "w", encoding="utf-8") as output_file:
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                output_file.write(line)

    # Open the input file in read mode and output file in write mode
    with open(inputfile, "r") as input_file, open(outputfile, "a") as output_file:
        # For each line, breakdow the SNPEff info field to retain only the effect entries
        # Iterate through the lines in the input file
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Get the EFF field (assuming it's always the eighth INFO field)
            info_field = fields[7]

            # Split the EFF field by commas to separete the SNPEff annotation
            aa, ac, af, snpeffs, *lo = info_field.split(';')

            # Check if any ReverseComplementedAlleles is present
            # in any SNPEff annotation entry and deal with it
            if (snpeffs == "ReverseComplementedAlleles"):
                _, snpeffann, = lo[0].split("EFF=")
            else:
                _, snpeffann, = snpeffs.split("EFF=")

            # Separete each SNPEff entry
            snpeffann_entries = snpeffann.split(',')

            # It only takes the first entry
            # Get the effect to use latter in keeping only relevant terms
            effect, _ = snpeffann_entries[0].split('(')

            # Re-assemble the EFF field with the simplified version
            # of the SNPEff annottion.
            fields[7] = f"{aa};{ac};{af};EFF={snpeffann_entries[0]}"

            # Re-assemble the entire line
            reassembled_line = '\t'.join(fields)

            if keeponlyterms:
                if any(x == effect for x in relevant_effect_terms):
                    output_file.write(reassembled_line + "\n")
            else:
                output_file.write(reassembled_line + "\n")

    return "file processed"


# simplify_snpeff_with_custom_annotation()
def simplify_snpeff_with_custom_annotation(
        file: str, custom_annotation: str,
        outfile: str = None, keeponlyterms: bool = False,
        method: str = "First"
        ) -> None:

    """
    This function simplifies SNPEff entries of a VCF file.
    This is method with custom annotations.
    For the moment, method don't do anything.
    """

    # Define input and output files
    inputfile, outputfile = parser_and_checker(file, outfile)

    # Check if a custom annotation file was provided
    if not custom_annotation:
        raise ValueError("custom_annotation must be provided")

    # This handles if keeponlyterms is set to True
    if keeponlyterms:
        # Define annotation effects terms we care about
        relevant_effect_terms = ['INTRON', 'SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING', 'INTERGENIC']

    # Copy the header from the inputfile to the ouputfile
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "w", encoding="utf-8") as output_file:
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                output_file.write(line)

    # Open the input file in read mode and output file in write mode
    with open(inputfile, "r") as input_file, open(outputfile, "a") as output_file:
        # For each line, breakdow the SNPEff info field to retain only the effect entries
        # Iterate through the lines in the input file
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Get the EFF field (assuming it's always the eighth INFO field)
            info_field = fields[7]

            # Split the EFF field by commas to separete the SNPEff annotation
            aa, ac, af, snpeffs, *lo = info_field.split(';')

            # Check if any ReverseComplementedAlleles is present
            # in any SNPEff annotation entry and deal with it
            if (snpeffs == "ReverseComplementedAlleles"):
                _, snpeffann, = lo[0].split("EFF=")
            else:
                _, snpeffann, = snpeffs.split("EFF=")

            # Separete each SNPEff entry
            snpeffann_entries = snpeffann.split(',')

            # It only takes the first entry
            # Get the effect to use latter in keeping only relevant terms
            effect, _ = snpeffann_entries[0].split('(')

            # Take SNPEff first entry and one CUSTOM if present
            custom_entries = [entry for entry in snpeffann_entries if entry.startswith(f"CUSTOM[{custom_annotation}]")]

            # Re-assemble the EFF field with the simplified version
            # of the SNPEff annottion.
            if (custom_entries == []):
                fields[7] = f"{aa};{ac};{af};EFF={snpeffann_entries[0]}"
            else:
                fields[7] = f"{aa};{ac};{af};EFF={snpeffann_entries[0]},{custom_entries[0]}"

            # Re-assemble the entire line
            reassembled_line = '\t'.join(fields)

            if keeponlyterms:
                if any(x == effect for x in relevant_effect_terms):
                    output_file.write(reassembled_line + "\n")
            else:
                output_file.write(reassembled_line + "\n")

    return "file processed"


def parseargs():
    parser = argparse.ArgumentParser("python simplify_snpeff.py", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="The the name of the file which the data are to be read from (annotated vcf from SNPEff)",
                        dest="file", required=True, type=str)
    parser.add_argument("-o", help="A character string naming a file",
                        dest="outfile", default=None, type=str)
    parser.add_argument("-f", help="Keep only relevant terms (it removes up(down)-stream SNPs)",
                        dest="keeponlyterms", default=False, type=bool)
    parser.add_argument("-m", help="Method to consolidate SNPEff annotation",
                        dest="method", default="First", type=str)
    parser.add_argument("-a", help="BED file with custom annotation",
                        dest="custom_annotation", default=None, type=str)
    return parser


def main(argv):
    """
    This is the main program definition.
    """
    parser = parseargs()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    file = args.file
    outfile = args.outfile
    keeponlyterms = args.keeponlyterms
    method = args.method
    custom_annotation = args.custom_annotation

    if custom_annotation is not None:
        result = simplify_snpeff_with_custom_annotation(
            file=file, outfile=outfile,
            keeponlyterms=keeponlyterms, method=method,
            custom_annotation=custom_annotation
        )
        print(result)

    else:
        result = simplify_snpeff_default(
            file=file, outfile=outfile,
            keeponlyterms=keeponlyterms, method=method
        )
        print(result)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        main(['-h'])
    else:
        main(sys.argv[1:])


# Majority-rule like effects filtering
# Create an empty dictionary to store the counts
# item_counts = {}

# # Loop through the list and count the occurrences of each item
# for item in snpeffann_effects_entries:
#     if item in item_counts:
#         item_counts[item] += 1
#     else:
#         item_counts[item] = 1

# # Use a threshold to define intron_variant to save
# # Calculate the total count of all items
# total_count = sum(item_counts.values())

# # Calculate the total count of the items you want to retain
# total_retain_count = sum(item_counts.get(item, 0) for item in items_to_retain)

# # Calculate the percentage for each item you want to retain
# item_percentages = {item: (item_counts.get(item, 0) / total_count) * 100 for item in items_to_retain}

# # Check if each item's percentage is >= 30%
# valid_items = [item for item in items_to_retain if item_percentages.get(item, 0) >= 30]

# # If all items meet the condition, retain them; otherwise, don't retain any
# if len(valid_items) > 0:
# #print(line)
# filtered_vcf_lines.append(line)

