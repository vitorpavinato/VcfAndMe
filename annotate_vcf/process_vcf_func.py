"""
process_vcf_func.py
This module contains a set of functions to process VCF files.
They are lower level functions, not for the user.
They are used in the highler level functions for different
instances of the implementations.py script.
"""

import re
from typing import Tuple, List


# Unpack vcf FIELDS items
def process_fields(
    fields: list
) -> Tuple[List[str], int, str, List[str]]:
    """
    Unpack vcf FIELDS items
    """
    # Unpack vcf FIELDS items
    chrom, pos, vcfid, ref, alt, qual, filtr, info = fields[:8]

    # Cast pos to integer to use with samtools faidx
    pos_int = int(pos)

    # Get genotypes
    genotypes = fields[9:]

    # Info list
    snp_fields_list = [chrom, pos, vcfid, ref, alt, qual, filtr]

    return snp_fields_list, pos_int, info, genotypes


# Unpack items in INFO
def process_info(
    info: str,
    sift4g_annotation: bool = False
) -> Tuple[List[str], str, str]:
    """
    Unpack items in INFO
    """

    # Unpack items in INFO
    aa_, ac_, af_, *annotations = info.strip().split(";")

    # Unpack subitems in INFO
    aa = aa_.split("=")[1]
    ac = ac_.split("=")[1]
    af = af_.split("=")[1]

    # Create info subitems list
    info_subitems_list = [aa, ac, af]

    snpeff = lof = nmd = sift = "NA"

    # Unpack items in annotations:
    for annotation in annotations:
        if re.search(r'EFF=.*', annotation):
            snpeff = annotation
        elif re.search(r'LOF=.*', annotation):
            lof = annotation
        elif re.search(r'NMD=.*', annotation):
            nmd = annotation
        elif re.search(r'SIFTINFO=.*', annotation):
            sift = annotation

    if sift4g_annotation:
        return info_subitems_list, snpeff, lof, nmd, sift
    # else:
    return info_subitems_list, snpeff, lof, nmd


# Count the number of alternative and reference genotypes
def count_genotypes(genotypes: list) -> Tuple[int, int]:
    """
    Count the number of alternative and reference genotypes
    in a list of genotypes and return a string with the counts
    """
    refcount = 0
    altcount = 0
    totalcount = 0

    for genotype in genotypes:
        if genotype != './.':
            if genotype == "0/0":
                refcount += 1
            elif genotype == "1/1":
                altcount += 1

    totalcount = refcount + altcount
    genotype_count_list = [refcount, altcount, totalcount]

    # Return a list with the counts
    return genotype_count_list


# Unpack items in info::snpeff (when present)
def get_snpeff_items(
    snpeff: str,
    custom_effect_name: str = None,
    new_custom_effect_name: str = None
) -> List[str]:
    """
    Create variables with NA to store snpeff items and avoid errors.
    Unpack items in snpeff.
    Return only needed items for the pipeline
    """
    refcodon, altcodon = "NA", "NA"
    maineffect = "NA"
    (
        effect_impact,
        functional_class,
        gene_name,
        transcript_biotype,
        gene_coding,
        transcript_id,
    ) = ("NA", "NA", "NA", "NA", "NA", "NA")

    # Unpack items in SNPeff info::eff
    _, effect_ = snpeff.strip().split("=")
    maineffect_, *customeffect_ = effect_.strip().split(",")

    # Unpack the main effect
    maineffect, eff = maineffect_.strip().split("(")
    eff, _ = eff.strip().split(")")

    # ALL SNPEff fields (gatk format)
    (
        effect_impact,
        functional_class,
        codon_change,
        aa_change,
        aa_length,
        gene_name,
        transcript_biotype,
        gene_coding,
        transcript_id,
        exon,
        genotype,
        *errors,
    ) = eff.strip().split("|")

    # This might fix empty values not replaced with NA
    if not bool(functional_class.strip()):
        functional_class = "NA"

    if not bool(transcript_biotype.strip()):
        transcript_biotype = "NA"

    # Unpack codons for functional classes
    if functional_class in ("SILENT", "MISSENSE", "NONSENSE"):
        refcodon, altcodon = codon_change.strip().split("/")

    # Unpack custom effects
    if custom_effect_name is not None:
        if len(customeffect_) > 0:
            ceff_, _ = customeffect_[0].strip().split("]")
            _, custom_effect_name_ = ceff_.strip().split("[")

            if new_custom_effect_name is not None:
                if new_custom_effect_name == custom_effect_name_:
                    customeffect = custom_effect_name_
                else:
                    customeffect = new_custom_effect_name
            else:
                customeffect = custom_effect_name_

            # Include the custom effect in the main effect
            maineffect = f"{maineffect}+{customeffect}"

    # Return only needed items for the pipeline
    snpeff_list = [refcodon, altcodon, maineffect, effect_impact, functional_class, gene_name, transcript_biotype, gene_coding, transcript_id]

    return snpeff_list


# Unpack items in info::sift (when present)
def get_sift4g_items(sift4g: str, threshold: float = 0.05) -> List[str]:
    """
    Create variables with NA to store sift4g items and avoid errors.
    Unpack items in sift4g.
    Return only needed items for the pipeline
    """

    # sift4g expect elements
    refaa, altaa = "NA", "NA"
    (
        transcript,
        geneid,
        genename,
        region,
        varianttype,
        siftscore,
        siftmedian,
        siftpred,
        deleteriousness,
    ) = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

    # Check if there is any SIFT4G info
    if sift4g != "NA":

        # Unpack items in info::eff
        _, siftinfo = sift4g.strip().split("=")

        # ALL SIFT fields
        (
            allele,
            transcript,
            geneid,
            genename,
            region,
            varianttype,
            aa,
            aaposition,
            siftscore,
            siftmedian,
            siftnumseqs,
            alleletype,
            siftpred,
        ) = siftinfo.strip().split("|")

        # Unpack items in sift4g aa item
        refaa, altaa = aa.strip().split("/")

        # Define deleteriousness status
        if siftscore != "NA":
            if float(siftscore) < threshold:
                deleteriousness = "deleterious"
            else:
                deleteriousness = "tolerated"

    # Return only needed items for the pipeline
    sift4g_list = [refaa, altaa, transcript, geneid, genename, region, varianttype, siftscore, siftmedian, siftpred, deleteriousness]

    return sift4g_list


# Unpack items in info::lof
def get_lof_items(lof: str) -> List[str]:
    """
    Create variables with NA to store lof items and avoid errors.
    Unpack items in lof.
    Return only needed items for the pipeline
    """

    # lof expect elements
    (
        genename,
        geneid,
        numbe_of_transcripts,
        perc_affected_transcripts
    ) = ("NA", "NA", "NA", "NA")

    # Check if there is any SIFT4G info
    if lof != "NA":

        # Unpack items in info::eff
        _, lof_ = lof.strip().split("=(")
        lof_, _ = lof_.strip().split(")")

        # ALL SIFT fields
        (
            genename,
            geneid,
            numbe_of_transcripts,
            perc_affected_transcripts
        ) = lof_.strip().split("|")

    # Return lof items
    lof_list = [genename, geneid, numbe_of_transcripts, perc_affected_transcripts]

    return lof_list


# Unpack items in info::nmd
def get_nmd_items(nmd: str) -> List[str]:
    """
    Create variables with NA to store lof items and avoid errors.
    Unpack items in nmd.
    Return only needed items for the pipeline
    """

    # lof expect elements
    (
        genename,
        geneid,
        numbe_of_transcripts,
        perc_affected_transcripts
    ) = ("NA", "NA", "NA", "NA")

    # Check if there is any SIFT4G info
    if nmd != "NA":

        # Unpack items in info::eff
        _, nmd_ = nmd.strip().split("=(")
        nmd_, _ = nmd_.strip().split(")")

        # ALL SIFT fields
        (
            genename,
            geneid,
            numbe_of_transcripts,
            perc_affected_transcripts
        ) = nmd_.strip().split("|")

    # Return lof items
    nmd_list = [genename, geneid, numbe_of_transcripts, perc_affected_transcripts]

    return nmd_list
