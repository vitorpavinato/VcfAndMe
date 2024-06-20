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
    sift4g_annotations: bool = False, sift_threshold: float = 0.05
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
        sift_threshold (float, optional): User defined version of the sift threshold. Defaults to 0.05.
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
            new_custom_effect_name=new_custom_effect_name,
            sift_threshold=sift_threshold
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


def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser("python vcf_to_tsv.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input vcf file name", dest="inputfile", required=True, type=str)
    parser.add_argument("-o", help="Output tsv file name", dest="outputfile", default=None, type=str)
    parser.add_argument("-r", help="Path to the reference genome of the vcf file", dest="reference", required=True, type=str)
    parser.add_argument("-s", help="Path to the samtools", dest="samtools_path", required=True, type=str)
    parser.add_argument("-f", help="Number of bases flanking each targeted SNP", dest="nflankinbps", default=3, type=int)
    parser.add_argument("-c", help="Custom effect name", dest="custom_effect_name", default=None, type=str)
    parser.add_argument("-n", help="New custom effect name", dest="new_custom_effect_name", default=None, type=str)
    parser.add_argument("-e", help="Input vcf with SIFT4G annotations", dest="sift4g_annotations", action="store_true")
    parser.add_argument("-d", help="User defined version of sift threshold for SIFT4G annotations", dest="sift_threshold", default=None, type=float)
    
    args = parser.parse_args()
    
    # Raise an error if -d is used without -e
    if args.sift_threshold is not None and not args.sift4g_annotations:
        parser.error("The -d flag requires the -e flag to be set. Please use -e when using -d.")
    
    # If -e is used without -d, set the default threshold
    if args.sift4g_annotations and args.sift_threshold is None:
        args.sift_threshold = 0.05
    
    return args


# Main program
def main(argv) -> None:
    """
    This is the main program definition.
    """

    args = parseargs()

    # Define input and output files
    inputfile = args.inputfile
    outputfile = args.outputfile
    reference = args.reference
    samtools_path = args.samtools_path
    nflankinbps = args.nflankinbps
    custom_effect_name = args.custom_effect_name
    new_custom_effect_name = args.new_custom_effect_name
    sift4g_annotations = args.sift4g_annotations
    sift_threshold = args.sift_threshold

    # Execute the function
    result = vcf_to_tsv(inputfile=inputfile, outputfile=outputfile,
                        reference=reference, samtools_path=samtools_path,
                        nflankinbps=nflankinbps, custom_effect_name=custom_effect_name,
                        new_custom_effect_name=new_custom_effect_name,
                        sift4g_annotations=sift4g_annotations,
                        sift_threshold=sift_threshold)

    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
