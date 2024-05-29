"""
process_vcf_func.py
This module contains a set of functions to process VCF files.
They are lower level functions, not for the user.
They are used in the highler level functions for different
instances of the implementations.py script.
"""


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
    aa_, ac_, af_, snpeff, *sift = info.strip().split(";")

    # Unpack subitems in INFO
    aa = aa_.split("=")[1]
    ac = ac_.split("=")[1]
    af = af_.split("=")[1]

    # Create info subitems list
    info_subitems_list = [aa, ac, af]

    if sift4g_annotation:
        return info_subitems_list, snpeff, sift
    # else:
    return info_subitems_list, snpeff


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

    # Unpack codons (for gatk format)
    # if effect == "SYNONYMOUS_CODING" or effect == "NON_SYNONYMOUS_CODING":
    if maineffect in ("SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING"):
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
    if len(sift4g) > 0:

        # Unpack items in info::eff
        _, siftinfo = sift4g[0].strip().split("=")
        siftinfo, *_ = siftinfo.strip().split(",")  # Bypass issue with tri-allelic SNPs

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
        if varianttype == "NONSYNONYMOUS":
            if siftscore == "NA":
                deleteriousness = "NA"
            else:
                if float(siftscore) < threshold:
                    deleteriousness = "deleterious"
                else:
                    deleteriousness = "tolerated"
        elif varianttype == "SYNONYMOUS":
            deleteriousness = "NA"

    # Return only needed items for the pipeline
    sift4g_list = [refaa, altaa, transcript, geneid, genename, region, varianttype, siftscore, siftmedian, siftpred, deleteriousness]

    return sift4g_list
