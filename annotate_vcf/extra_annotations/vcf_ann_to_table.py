"""
Simple program to convert SNPEff interval annotation only vcf to a table.
This program a supplementary module of VcfAndMe.
It can be used to process any extra annotation made with SNPEff
using (any) interval (it reads .BED) file.
The extra annotations as tsv file can be easily integrated
into an existing tsv file using Pandas, for exampe.
"""

import re
import sys
import argparse
import warnings
from utils import *


def parse_vcf_ann(
    inputfile: str, outputfile: str, annotation_name: str = None
) -> None:
    """
    Parse SNPEff vcf annotation to a table.
    """

    # Input file check
    inputfile = check_input_file(inputfile)

    # Output file check
    outputfile = output_file_name(inputfile, outputfile)

    # Check if the annotation name is provided
    if annotation_name == "" or annotation_name is None:
        warnings.warn("Annotation name is empty or None. Using the annotation name provided by SNPEff", category=UserWarning)

    # Initialize the header
    header = ["chrom", "position", "custom_annotation"]
    header = "\t".join(str(item) for item in header)

    # Initialize the a list to store the vcf lines
    vcflines = []

    with open(inputfile, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Split the line into fields using tab as the delimiter
            vcffields = line.strip().split('\t')

            # Get the INFO field (assuming it's always the eighth INFO field)
            info_field = vcffields[7]

            # Split the INFO field by commas to separete the SNPEff annotation
            aa, ac, af, *annotations = info_field.split(';')

            # Initialize the custom effect
            # Since it is a custom only SNP annotation, it is ok to
            # set the customeffect to NA
            customeffect = 'NA'

            # Iterate over the annotations list and check for patterns
            if len(annotations) != 0:
                for annotation in annotations:
                    if re.search(r'EFF=.*', annotation):
                        # Separate the annotation from EFF
                        _, snpeffentries = annotation.split("EFF=")

                        # Take the first custom effect
                        # Consistent with simplify_snpeff.py behavior
                        snpeffentry = snpeffentries.split(',')[0]

                        # Unpack the custom effect
                        _, customeffect_ = snpeffentry.strip().split("[")
                        customeffect, _ = customeffect_.strip().split("]")

                    elif re.search(r'ReverseComplementedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        # vcffields[2] = vcffields[2] + "+" + annotation
                        pass

                    elif re.search(r'SwappedAlleles', annotation):
                        # Take any of the ReverseComplementedAlleles and
                        # add it to the vcf ID fields
                        # vcffields[2] = vcffields[2] + "+" + annotation
                        pass

            vcfline = [vcffields[0], vcffields[1], customeffect]
            vcflines.append(vcfline)

    if annotation_name is not None:
        for i in range(len(vcflines)):
            if vcflines[i][2] != 'NA':
                vcflines[i][2] = annotation_name

    # Write the vcf lines to the output file:
    write_tsv_file(vcflines, header, outputfile)

    return "Terminated successfully"


# Define the command line interface
def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser("python vcf_ann_to_table.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input vcf file name", dest="inputfile", required=True, type=str)
    parser.add_argument("-o", help="Output tsv file name", dest="outputfile", default=None, type=str)
    parser.add_argument("-n", help="New custom annotation name", dest="annotation_name", default=None, type=str)
    return parser


# Main function
def main(argv) -> None:
    """
    This is the main program definition.
    """
    parser = parseargs()
    if argv[-1] == "":
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define input and output files
    inputfile = args.inputfile
    outputfile = args.outputfile
    annotation_name = args.annotation_name

    # Execute the function
    result = parse_vcf_ann(inputfile=inputfile,
                           outputfile=outputfile,
                           annotation_name=annotation_name)

    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
