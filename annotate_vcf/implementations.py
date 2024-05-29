"""
This contains different implementations of the vcf_to_tsv main functions.
"""

from typing import Tuple, List
from process_vcf_func import *
from mutational_context_func import get_mutational_context


# Process SNPEff and SIFT4G annotations
def processes_snpeff_sift4g_vcf(
    inputfile: str,
    reference: str, samtools: str, nflankinbps: int,
    custom_effect_name: str = None, new_custom_effect_name: str = None
) -> Tuple[list, list, list]:
    """
    Railroad pattern #1: The vcf includes annotations from SNPEff and SIFT4G
    """

    print("Start processing the SNPEff and SIFT4G annotated VCF file. It might take a while...")
    with open(inputfile, "r", encoding="utf-8") as input_file:

        # Initialize variables used internally
        current_block = []
        lines_in_block = []
        previous_position = 0
        line_trck = 0

        # Initialize variables that will be returned
        list_block = []
        list_lines_in_block = []
        vcf_lines = []

        # Read each line in the vcf file and save to a vcf_lines
        for line in input_file:
            if line.startswith("##") or line.startswith("#"):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            # Check if the line has at least 8 elements
            if len(fields) < 8:
                continue

            # Unpack vcf FIELDS items
            snp_fields_list, pos_int, info, genotypes = process_fields(fields=fields)

            # Unpack vcf INFO items
            info_subitems_list, snpeff_, sift_ = process_info(info=info, sift4g_annotation=True)

            # Count genotypes
            genotype_count_list = count_genotypes(genotypes)

            # Get the mutational context
            mutational_context_list = get_mutational_context(
                chrom=snp_fields_list[0],
                pos=pos_int,
                refallele=snp_fields_list[3],
                altallele=snp_fields_list[4],
                reference=reference,
                samtools=samtools,
                nflankinbps=nflankinbps
            )

            # Unpack snpeff INFO items
            snpeff_list = get_snpeff_items(
                snpeff=snpeff_,
                custom_effect_name=custom_effect_name,
                new_custom_effect_name=new_custom_effect_name
            )

            # Unpack SIFT4G INFO items
            sift4g_list = get_sift4g_items(sift_, threshold=0.05)

            # This part is for fixing basis in flaking bases string
            # if SNPs are close (less than nflankinbps)
            # This keeps track of block of positions with near SNPs
            current_position = pos_int

            if (
                current_position - previous_position
            ) < nflankinbps:  # (maybe nflankinbps - 1)
                if not current_block:
                    current_block.append(previous_position)
                    lines_in_block.append(line_trck - 1)
                current_block.append(current_position)
                lines_in_block.append(line_trck)
            else:
                current_block = []
                lines_in_block = []

            if len(current_block) > 1:
                if current_block not in list_block:
                    list_block.append(current_block)
                if lines_in_block not in list_lines_in_block:
                    list_lines_in_block.append(lines_in_block)

            previous_position = current_position

            # Create a re-usable list to store the needed information
            line_list = snp_fields_list + info_subitems_list + genotype_count_list + mutational_context_list + snpeff_list + sift4g_list

            # Add the line to the list of lines
            vcf_lines.append(line_list)
            line_trck += 1

    return vcf_lines, list_lines_in_block, list_block


# Process SNPEff ONLY annotations
def processes_snpeff_vcf(
    inputfile: str,
    reference: str, samtools: str, nflankinbps: int,
    custom_effect_name: str = None, new_custom_effect_name: str = None
) -> Tuple[list, list, list]:
    """
    Railroad pattern #2: The vcf includes annotations from SNPEff only
    """

    print("Start processing the SNPEff only annotated VCF file. It might take a while...")
    with open(inputfile, "r", encoding="utf-8") as input_file:

        # Initialize variables used internally
        current_block = []
        lines_in_block = []
        previous_position = 0
        line_trck = 0

        # Initialize variables that will be returned
        list_block = []
        list_lines_in_block = []
        vcf_lines = []

        # Read each line in the vcf file and save to a vcf_lines
        for line in input_file:
            if line.startswith("##") or line.startswith("#"):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            # Check if the line has at least 8 elements
            if len(fields) < 8:
                continue

            # Unpack vcf FIELDS items
            snp_fields_list, pos_int, info, genotypes = process_fields(fields=fields)

            # Unpack vcf INFO items
            info_subitems_list, snpeff_ = process_info(info=info, sift4g_annotation=False)

            # Count genotypes
            genotype_count_list = count_genotypes(genotypes)

            # Get the mutational context
            mutational_context_list = get_mutational_context(
                chrom=snp_fields_list[0],
                pos=pos_int,
                refallele=snp_fields_list[3],
                altallele=snp_fields_list[4],
                reference=reference,
                samtools=samtools,
                nflankinbps=nflankinbps
            )

            # Unpack snpeff INFO items
            snpeff_list = get_snpeff_items(
                snpeff=snpeff_,
                custom_effect_name=custom_effect_name,
                new_custom_effect_name=new_custom_effect_name
            )

            # This part is for fixing basis in flaking bases string
            # if SNPs are close (less than nflankinbps)
            # This keeps track of block of positions with near SNPs
            current_position = pos_int

            if (
                current_position - previous_position
            ) < nflankinbps:  # (maybe nflankinbps - 1)
                if not current_block:
                    current_block.append(previous_position)
                    lines_in_block.append(line_trck - 1)
                current_block.append(current_position)
                lines_in_block.append(line_trck)
            else:
                current_block = []
                lines_in_block = []

            if len(current_block) > 1:
                if current_block not in list_block:
                    list_block.append(current_block)
                if lines_in_block not in list_lines_in_block:
                    list_lines_in_block.append(lines_in_block)

            previous_position = current_position

            # Create a re-usable list to store the needed information
            line_list = snp_fields_list + info_subitems_list + genotype_count_list + mutational_context_list + snpeff_list

            # Add the line to the list of lines
            vcf_lines.append(line_list)
            line_trck += 1

    return vcf_lines, list_lines_in_block, list_block


def snpeff_sift4g_header() -> str:
    """
    This function returns the header of the output file
    """
    header = [
        "chrom", "pos", "aainfo", "ref", "alt", "qual", "filtr", "aa", "ac", "af",
        "refcount", "altcount", "totalcount", "refcontext", "altcontext",
        "refcontext_complrev", "altcontext_complrev",
        "refcodon", "altcodon", "maineffect", "snpeff_effimpact", "snpeff_funclass",
        "snpeff_genename", "snpeff_trnscbiotype", "snpeff_genecoding", "snpeff_trnscid",
        "refaa", "altaa", "sift_trnscid", "sift_geneid", "sift_genename", "sift_region",
        "sift_vartype", "sifts_core", "sift_median", "sift_pred", "deleteriousness"
    ]
    header = "\t".join(str(item) for item in header)

    return header


def snpeff_header() -> str:
    """
    This function returns the header of the output file
    """
    header = [
        "chrom", "pos", "aainfo", "ref", "alt", "qual", "filtr", "aa", "ac", "af",
        "refcount", "altcount", "totalcount", "refcontext", "altcontext",
        "refcontext_complrev", "altcontext_complrev",
        "refcodon", "altcodon", "maineffect", "snpeff_effimpact", "snpeff_funclass",
        "snpeff_genename", "snpeff_trnscbiotype", "snpeff_genecoding", "snpeff_trnscid"
    ]
    header = "\t".join(str(item) for item in header)

    return header
