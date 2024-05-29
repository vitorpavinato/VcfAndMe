"""
Implementation of the annotation pipeline found on Jupyter Notebooks.
"""

import os
import argparse
import sys
from pathlib import Path
from annotation_class import SIFT4GAnnotation, SNPEFFAnnotation


def parseargs():
    parser = argparse.ArgumentParser("python pipeline_to_annotate_vcf.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="input vcf file not annotated",
                        dest="input_file", required=True, type=str)
    parser.add_argument("-d", help="SNPEff database",
                        dest="snpeff_database", required=True, type=str)
    parser.add_argument("-b", help="SNPEff interval BED file",
                        dest="snpeff_interval_file", default=None, type=str)
    parser.add_argument("-o", help="output folder for SNPEff annotations",
                        dest="snpeff_output_folder", required=True, type=str)
    parser.add_argument("-s", help="SIFT4G database",
                        dest="sift4g_database", default=None, type=str)
    parser.add_argument("-f", help="Output folder for SIFT4G annotations",
                        dest="sift4g_output_folder", default=None, type=str)
    return parser


def main(argv) -> None:
    """
    Main function
    """
    parser = parseargs()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    INPUT_FILE = args.input_file
    SNPEFF_DATABASE = args.snpeff_database
    SNPEFF_INTERVAL_FILE = args.snpeff_interval_file
    SNPEFF_OUTPUT_FOLDER = args.snpeff_output_folder
    SIFT4G_DATABASE = args.sift4g_database
    SIFT4G_OUTPUT_FOLDER = args.sift4g_output_folder

    # Check if the input file exists
    if not Path(INPUT_FILE).exists():
        raise ValueError("Input file does not exist")

    # Run the SNPEFF annotation
    # This is the minimal behavior of this pipeline

    # Create the SNPEFFAnnotation object
    snpeff_annotation = SNPEFFAnnotation(
        input_file=INPUT_FILE,
        database=SNPEFF_DATABASE,
        output_folder=SNPEFF_OUTPUT_FOLDER
    )

    # Check if the interval file was provided
    if SNPEFF_INTERVAL_FILE not in [None, ""]:

        # Get the path and filename
        path, filename = os.path.split(SNPEFF_INTERVAL_FILE)

        # Define the basename for the outputs from the filename
        custom_annotation, _ = filename.strip().split('.bed')

        # Run the SNPEFF annotation with custom annotations
        snpeff_annotation.annotate_with_intervals(interval_file=SNPEFF_INTERVAL_FILE)

        # Run simplify_snpeff annotation with custom annotations
        snpeff_annotation.simplify_snpeff_annotations(custom_annotation=custom_annotation)

    else:
        # Run the SNPEFF annotation without custom annotations
        snpeff_annotation.annotate()

        # Run simplify_snpeff annotation without custom annotations
        snpeff_annotation.simplify_snpeff_annotations()

    # Run the optional SIFT4G annotation
    if SIFT4G_DATABASE not in [None, ""]:
        if SIFT4G_OUTPUT_FOLDER is None:
            raise ValueError("sift4g output folder must be provided")

        # Create the SIFT4GAnnotation object
        sift4g_annotation = SIFT4GAnnotation(
            input_file=snpeff_annotation.output_file_simplified,
            database=SIFT4G_DATABASE,
            output_folder=SIFT4G_OUTPUT_FOLDER
        )

        # Run SIFT4G annotation
        sift4g_annotation.annotate()


if __name__ == "__main__":

    if len(sys.argv) < 2:
        main(['-h'])
    else:
        main(sys.argv[1:])
