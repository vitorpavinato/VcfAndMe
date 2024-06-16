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
import re


# Function to check if a file exists
def parser_and_checker(file: str, outfile: str = None) -> bool:
    """
    Simple function to check:
    - if the input file was provided by the user and if it exists;
    - if the input file has something other than the header; 
    - if the input file has at least 4 elements in the INFO fields;
    - the 4 elements in the INFO field should be::
        - AA, AC, AF, and EFF
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

    # Process the first input line
    first_line = None

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

        # Check if the INFO fields:
        # Should have at least 4 elements in the list: AC, AF, AA, EFF
        if len(info_field) >= 4:
            if not (any(re.search(r'AA=.*', element) for element in info_field) and
                    any(re.search(r'AC=.*', element) for element in info_field) and
                    any(re.search(r'AF=.*', element) for element in info_field) and
                    any(re.search(r'EFF=.*', element) for element in info_field)):
                raise ValueError("Input is not supported by this script! It should have at least AA, AC, AF,and EFF")

        else:
            # Handle the case when the INFO field has not enough elements
            raise ValueError(
                "Input file has not enough information in the INFO field...")

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
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "a", encoding="utf-8") as output_file:
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
            aa, ac, af, *annotations = info_field.split(';')

            # Start re-assembling INFO field
            fields[7] = f"{aa};{ac};{af}"

            # Iterate over the annotations list and check for patterns
            # This handle expected elements in annotation more explicitly
            if len(annotations) != 0:
                for annotation in annotations:
                    if re.search(r'EFF=.*', annotation):
                        # Separate the annotation for EFF
                        _, effentries = annotation.split("EFF=")

                        # Take the first EFF effect
                        effentry = effentries.split(',')[0]

                        # Add the simplified EFF field to the INFO field
                        fields[7] = f"{fields[7]};EFF={effentry}"

                    elif re.search(r'LOF', annotation):
                        # Separate the annotation for LOF
                        _, lofentries = annotation.split("LOF=")

                        # Take the first LOF effect
                        # Consistent with EFF effect above
                        lofentry = lofentries.split(',')[0]

                        # Add the simplified LOF field to the INFO field
                        fields[7] = f"{fields[7]};LOF={lofentry}"

                    elif re.search(r'NMD', annotation):
                        # Separate the annotation for NMD
                        _, nmdentries = annotation.split("NMD=")

                        # Take the first NMD effect
                        # Consistent with EFF effect above
                        nmdentry = nmdentries.split(',')[0]

                        # Add the simplified NMD field to the INFO field
                        fields[7] = f"{fields[7]};NMD={nmdentry}"

                    elif re.search(r'ReverseComplementedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        fields[2] = fields[2] + "+" + annotation
                        # pass

                    elif re.search(r'SwappedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        fields[2] = fields[2] + "+" + annotation
                        # pass

            # It only takes the first entry
            # Get the effect to use latter in keeping only relevant terms
            effect, _ = effentry.split('(')

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
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "a", encoding="utf-8") as output_file:
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
            aa, ac, af, *annotations = info_field.split(';')

            # Start re-assembling INFO field
            fields[7] = f"{aa};{ac};{af}"

            # Iterate over the annotations list and check for patterns
            # This handle expected elements in annotation more explicitly
            if len(annotations) != 0:
                for annotation in annotations:
                    if re.search(r'EFF=.*', annotation):
                        # Separate the annotation for EFF
                        _, effentries = annotation.split("EFF=")

                        # Split the EFF entries by commas
                        effentries = effentries.split(',')

                        # Take the first EFF effect to simplify it
                        fields[7] = f"{fields[7]};EFF={effentries[0]}"

                    elif re.search(r'LOF', annotation):
                        # Separate the annotation for LOF
                        _, lofentries = annotation.split("LOF=")

                        # Take the first LOF effect
                        # Consistent with EFF effect above
                        lofentry = lofentries.split(',')[0]

                        # Add the simplified LOF field to the INFO field
                        fields[7] = f"{fields[7]};LOF={lofentry}"

                    elif re.search(r'NMD', annotation):
                        # Separate the annotation for NMD
                        _, nmdentries = annotation.split("NMD=")

                        # Take the first NMD effect
                        # Consistent with EFF effect above
                        nmdentry = nmdentries.split(',')[0]

                        # Add the simplified NMD field to the INFO field
                        fields[7] = f"{fields[7]};NMD={nmdentry}"

                    elif re.search(r'ReverseComplementedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        fields[2] = fields[2] + "+" + annotation
                        # pass

                    elif re.search(r'SwappedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        fields[2] = fields[2] + "+" + annotation
                        # pass

            # It only takes the first entry
            # Get the effect to use latter in keeping only relevant terms
            effect, _ = effentries[0].split('(')

            # Take SNPEff first entry and one CUSTOM if present
            custom_entries = [entry for entry in effentries if entry.startswith(f"CUSTOM[{custom_annotation}]")]

            # Re-assemble the EFF with the CUSTOM field
            if custom_entries != []:
                # Add the simplified CUSTOM field to the INFO field
                fields[7] = f"{fields[7]},{custom_entries[0]}"

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
