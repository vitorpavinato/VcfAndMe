"""
Program converts an annotated and simplified vcf file to a .tsv table
"""

import argparse
import sys
from vcf_to_tsv_func import *
from implementations import *
from mutational_context_func import fix_mutational_context


# vcf_to_tsv.py railroad pattern implementation
def vcf_to_tsv(
    inputfile: str, outputfile: str,
    reference: str, samtools_path: str, nflankinbps: int,
    custom_effect_name: str = None, new_custom_effect_name: str = None,
    sift4g_annotations: bool = False
) -> None:
    """vcf_to_tsv.py railroad pattern implementation.
    This function takes a vcf file and converts it to a .tsv table.
    The vcf might have SNPEFF annotations or SIFT4G annotations.

    Args:
        inputfile (str): Input vcf file name
        outputfile (str): Output tsv file name
        reference (str): Path to the reference genome of the vcf file
        samtools_path (str): Path to the samtools
        nflankinbps (int): Number of bases flanking each targeted SNP
        custom_effect_name (str, optional): Custom effect name. Defaults to None.
        new_custom_effect_name (str, optional): New custom effect name. Defaults to None.
        sift4g_annotations (bool, optional): Input vcf with SIFT4G annotations. Defaults to False.
    """

    # Input file check
    inputfile = check_input_file(inputfile)

    # Output file check
    outputfile = output_file_name(inputfile, outputfile)

    # Samtools path check
    samtools = check_samtools_path(samtools_path)

    # Reference file check
    reference = check_reference_genome_file(reference)

    # Create .fai file if it doesn't exist
    create_faidx(reference, samtools)

    # Here is the railroad pattern implementation
    if sift4g_annotations:
        vcf_lines, list_lines_in_block, list_block = processes_snpeff_sift4g_vcf(
            inputfile=inputfile, reference=reference,
            samtools=samtools, nflankinbps=nflankinbps,
            custom_effect_name=custom_effect_name,
            new_custom_effect_name=new_custom_effect_name
        )
        header = snpeff_sift4g_header()

    else:
        vcf_lines, list_lines_in_block, list_block = processes_snpeff_vcf(
            inputfile=inputfile, reference=reference,
            samtools=samtools, nflankinbps=nflankinbps,
            custom_effect_name=custom_effect_name,
            new_custom_effect_name=new_custom_effect_name
        )
        header = snpeff_header()

    # Execute the rest of the function.
    vcf_lines = fix_mutational_context(
        list_lines_in_block=list_lines_in_block,
        list_block=list_block,
        vcf_lines=vcf_lines,
        nflankinbps=nflankinbps
    )

    # Write .tsv file
    write_tsv_file(vcf_lines, header, outputfile)

    return "vcf to tsv processed"


# Define the command line arguments
def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser("python vcf_to_table.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input vcf file name", dest="inputfile", required=True, type=str)
    parser.add_argument("-o", help="Output tsv file name", dest="outputfile", default=None, type=str)
    parser.add_argument("-r", help="Path to the reference genome of the vcf file", dest="reference", required=True, type=str)
    parser.add_argument("-s", help="Path to the samtools", dest="samtools_path", required=True, type=str)
    parser.add_argument("-f", help="Number of bases flanking each targeted SNP", dest="nflankinbps", default=3, type=int)
    parser.add_argument("-c", help="Custom effect name", dest="custom_effect_name", default=None, type=str)
    parser.add_argument("-n", help="New custom effect name", dest="new_custom_effect_name", default=None, type=str)
    parser.add_argument("-e", help="Input vcf with SIFT4G annotations", dest="sift4g_annotations", default=False, action="store_true")
    return parser


# Main program
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
    reference = args.reference
    samtools_path = args.samtools_path
    nflankinbps = args.nflankinbps
    custom_effect_name = args.custom_effect_name
    new_custom_effect_name = args.new_custom_effect_name
    sift4g_annotations = args.sift4g_annotations

    # Execute the function
    result = vcf_to_tsv(inputfile=inputfile, outputfile=outputfile,
                        reference=reference, samtools_path=samtools_path,
                        nflankinbps=nflankinbps, custom_effect_name=custom_effect_name,
                        new_custom_effect_name=new_custom_effect_name,
                        sift4g_annotations=sift4g_annotations)

    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
