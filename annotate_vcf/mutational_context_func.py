"""
This module contains the function to fix mutational context
of processed vcf lines.
"""


import subprocess
from typing import Tuple, List


# Get reverse complement of a DNA sequence
def get_reversed_complementary_strand(sequence: str) -> str:
    """
    Find the complementary strand and reverse to 5' to 3'.
    """
    # This is a dictionary with the complement nucleotides
    nt_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}

    # Cast the input sequence to a list
    sequence = list(sequence)

    # Obtain the complementary strand
    complementary = [nt_pairs[nt] for nt in sequence]

    # Reverse the complementary strand
    reverse = list(complementary[::-1])

    # Return the reversed complementary strand
    return "".join(reverse)


# Get the mutational context
def get_mutational_context(
    chrom: str,
    pos: int,
    refallele: str,
    altallele: str,
    reference: str,
    samtools: str,
    nflankinbps: int = 3
) -> Tuple[str, str]:
    """
    Get the mutational context.
    """

    # Initialize the mutational context to avoid errors.
    refcontext, altcontext = "NA", "NA"

    # Get bases before the SNP
    nt_bf = subprocess.run(
        [
            samtools,
            "faidx",
            reference,
            (chrom + ":" + str(pos - nflankinbps) + "-" + str(pos - 1)),
        ],
        capture_output=True, check=False
    )

    # Get bases after the SNP
    nt_af = subprocess.run(
        [
            samtools,
            "faidx",
            reference,
            (chrom + ":" + str(pos + 1) + "-" + str(pos + nflankinbps)),
        ],
        capture_output=True, check=False
    )

    # Remove new line characters
    _, flkng_bf, *_ = str(nt_bf.stdout).strip().split("\\n")
    _, flkng_af, *_ = str(nt_af.stdout).strip().split("\\n")

    # Get the mutational context
    refcontext = flkng_bf.upper() + refallele + flkng_af.upper()
    altcontext = flkng_bf.upper() + altallele + flkng_af.upper()

    # Get the complementary strand for each allele context
    refcontext_complrev = get_reversed_complementary_strand(refcontext)
    altcontext_complrev = get_reversed_complementary_strand(altcontext)

    # Return the mutational context as a list
    # In the same order as in the header function.
    mutational_context_list = [refcontext, altcontext,
                               refcontext_complrev, altcontext_complrev]

    # Return the mutational context
    return mutational_context_list


def fix_mutational_context(
        list_lines_in_block: List[List[int]], 
        list_block: List[int], 
        vcf_lines: List[List[str]], 
        nflankinbps: int
) -> List[List[str]]:
    """
    This function fix the mutational context of processed vcf lines.
    It runs at the end, when all lines were processed. It fixes the
    refcontext and altcontext of closest SNPs.
    """

    print("Fixing ref and alt mutations on flanking bases. It might take a while...")
    # for block_lines, block_pos in zip(new_list_lines_in_block, new_list_block):
    for block_lines, block_pos in zip(list_lines_in_block, list_block):
        blocksize = len(block_lines)
        if blocksize > 1:
            for i, current_element in enumerate(block_lines):
                elements_before = block_lines[:i]
                elements_after = block_lines[i + 1:]
                current_pos = block_pos[i]
                pos_before = block_pos[:i]
                pos_after = block_pos[i + 1:]

                # Get the flanking bases that need to be updated
                current_refcontext = list(vcf_lines[current_element][13])  # refcontext
                current_altcontext = list(vcf_lines[current_element][14])  # altcontext

                # Check if there are multiple elements before or after
                if len(elements_before) >= 1:
                    for j, eb in enumerate(elements_before):
                        dist = abs(current_pos - pos_before[j])
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_before = vcf_lines[eb][3]  # ref before
                            alt_before = vcf_lines[eb][4]  # alt before
                            current_refcontext[nflankinbps - dist] = ref_before
                            current_altcontext[nflankinbps - dist] = alt_before

                if len(elements_after) >= 1:
                    for j, ea in enumerate(elements_after):
                        dist = abs(current_pos - pos_after[j])
                        print(dist)
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_after = vcf_lines[ea][3]  # ref after
                            alt_after = vcf_lines[ea][4]  # alt after
                            current_refcontext[nflankinbps + dist] = ref_after
                            current_altcontext[nflankinbps + dist] = alt_after

                # Here update the flanking sequences in the corresponding line
                vcf_lines[current_element][13] = "".join(current_refcontext)
                vcf_lines[current_element][14] = "".join(current_altcontext)

                # Here update the complementary flanking sequences
                vcf_lines[current_element][15] = get_reversed_complementary_strand(vcf_lines[current_element][13])
                vcf_lines[current_element][16] = get_reversed_complementary_strand(vcf_lines[current_element][14])

    print("Mutational context fixed!")
    return vcf_lines
