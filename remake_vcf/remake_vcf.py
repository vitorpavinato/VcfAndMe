'''
Remake a snp-site vcf to a standard format.
I tried to implement a linux principle of having
a tool to do a well a specific job. This tool is a
bit more comphrehensive as it do some housekeeping 
work to make a proper vcf file.
'''
import sys
import argparse
import subprocess  # this is being used to run samtools
from pathlib import Path  # This is being used to check if the file exists
from typing import List


# Check if .fai exists and create it with samtools if it doesn't
def create_samtools_faix(reference_file: Path, samtools_path: Path) -> None:
    """
    Check the existence of a reference file faidx
    or create it with samtools
    """

    # Check if .fai associated to a reference genome file in fasta exists    
    samtools = samtools_path
    reference_fai_file = reference_file.with_suffix(reference_file.suffix + ".fai")
    if not reference_fai_file.exists():
        make_faidx = subprocess.run([samtools, "faidx", reference_file],
                                    capture_output=True, check=True)


def remake_vcf_header(
    input_file: Path, reference_file: Path, chrom_name: str,
    chrom_length: int, output_file: Path
) -> None:
    """
    Remake a vcf header to a standard format compatible with GATK
    """
    # Process header lines
    header_lines = []
    with open(input_file, 'r', encoding='utf-8') as input_file_obj:
        # Process lines starting with "##" or "#"
        for line in input_file_obj:
            if line.startswith("##") or line.startswith("#"):
                header_lines.append(line)

    # Add new lines and modify existing lines
    modified_header_lines = header_lines.copy()  # Create a copy of the original header lines

    # Remove two lines but add then later in the file processing
    # That is the reason to save then as variable for later
    # at same time I am using pop() method.
    chrm_line_info = modified_header_lines.pop(1)
    format_str_inf = modified_header_lines.pop(1)
    colum_names_line = modified_header_lines.pop(1)

    # Add new line
    modified_header_lines.append("##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n")

    # Add FORMAT line removed above
    modified_header_lines.append(format_str_inf)

    # Add multiple new lines
    new_lines = [
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n",
        #"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
    ]
    modified_header_lines.extend(new_lines)

    # Add a new and modified contig/chromosome line information
    modified_header_lines.append(f"##contig=<ID={chrom_name},length={chrom_length}>\n")

    # Add line with reference value
    modified_header_lines.append(f"##reference=file://{reference_file}\n")

    # Add column names line
    modified_header_lines.append(colum_names_line)

    # Write the modified header lines to the output file or perform further manipulations
    # For example, writing to the output file:
    with open(output_file, 'w', encoding='utf-8') as output_file_obj:
        for line in modified_header_lines:
            output_file_obj.write(line)


def remake_vcf(
    input_file: Path, reference_file: Path = None,
    chrom_name: str = None, chrom_length: int = None,
    output_file: Path = None, samtools_path: Path = None
) -> None:
    """
    Remake a vcf to a standard format compatible with GATK
    It assumes that the vcf file to be remade is a snp-site
    vcf for one chromosome. Chromosomes will be named as chr2L,
    chr2R, etc., since UCSC references have chromosomes prefix
    name as chr. In this case, the reference file and the
    chrom_name must be provided for the same chromosome as
    in the vcf file.
    """

    # Check if all required arguments are provided
    # This will raise an error if any of the required arguments are None

    # This is for the reference file
    if reference_file is None:
        raise ValueError("reference_file must be provided")

    # This is for the chrom_name
    if chrom_name is None:
        raise ValueError("chrom_name must be provided")

    # This is for the output_file
    if output_file is None:
        raise ValueError("output_file must be provided")

    # This is for the samtools_path
    if samtools_path is None:
        raise ValueError("samtools_path must be provided")

    # This is for the chrom_length
    if chrom_length is None:
        chrom_length = 9999999999

    # Check the existence of a reference file faidx
    create_samtools_faix(reference_file, samtools_path)
    print("Faidx created...")

    # Run remake_vcf_header
    remake_vcf_header(input_file, reference_file,
                      chrom_name, chrom_length,
                      output_file)
    print("New header created...")

    # Process remaining lines
    with open(input_file, "r", encoding="utf-8") as input_file_obj, open(output_file, "a", encoding="utf-8") as output_file_obj:
        for line in input_file_obj:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            # Skip empty lines or lines with less than 9 fields
            if len(fields) < 9:
                continue

            # Replace or add chromosome name
            # This should take care of chromosome names like chr2L, chr2R, chrX, etc
            # If a wrong name is provided, it will be replaced.
            if fields[0].startswith("chr"):
                if fields[0] != chrom_name:
                    print("Chromsome name is not the same as provided. It is going to be replaced.")
                    fields[0] = chrom_name

            # This should take care of chromosome names like 2L, 2R, X, etc
            # If a wrong chromosome was provided, it will be replaced.
            elif fields[0].endswith(("L", "R")) or fields[0] == "X":
                if fields[0] != chrom_name.lstrip("chr"):
                    print("Chromsome name is not the same as provided. It is going to be replaced.")
                    fields[0] = chrom_name

            # And this should take care of unnamed chromosomes with "." or 1
            # This should be from the snp-site vcf
            elif fields[0] in (".", "1"):
                print(f"Chromosome {chrom_name} will be added to the SNP.")
                fields[0] = chrom_name

            # Get vcf ref and alt bases
            vcf_candidate_ref_base = fields[3]
            vcf_candidate_alt_bases = fields[4].split(",")

            # Create a list with all bases: reference and alternatives
            # The assigned snps-sites reference should the the first
            vcf_candidate_bases = [vcf_candidate_ref_base] + vcf_candidate_alt_bases

            # Check reference base with samtools faidx
            pos = int(fields[1])  # First, cast pos to integer

            # Call samtools faidx to check reference base
            genome_ref_base = subprocess.run([
                samtools_path,
                "faidx",
                reference_file,
                (chrom_name + ":" + str(pos) + "-" + str(pos)),
            ],
                                    capture_output=True,
                                    check=False)

            # Get the reference base as a variable
            _, genome_ref_base, *_ = str(genome_ref_base.stdout).strip().split("\\n")
            
            # This will make sure that masked genome will work
            genome_ref_base = genome_ref_base.upper()

            # Check if the genome reference base is present in the list of vcf candidate bases
            # This might deal with assertiment bias on the reference genome
            if genome_ref_base not in vcf_candidate_bases:
                print("The reference base is not present in the list of vcf candidate bases. It will be appended.")
                vcf_candidate_bases.append(genome_ref_base)

            # Create a list of indexes based on the original order of vcf_candidate_bases elements
            vcf_candidate_bases_idx = [i for i, x in enumerate(vcf_candidate_bases)]
            
            # Get the index of the missing data character '*'
            idx_of_missing_data = vcf_candidate_bases.index("*")

            # Convert the missing data code to the standard format
            fields[9:] = ["./." if x == str(idx_of_missing_data) else x for x in fields[9:]]

            # Then remove the missing data character index from vcf_candidate_bases_idx
            vcf_candidate_bases_idx.remove(idx_of_missing_data)

            # Get the index of the genome reference base on the list of vcf candidate bases
            idx_of_ref_base = vcf_candidate_bases.index(genome_ref_base)

            # Create a list of re-ordered indexes where the first value is the index of the reference base
            vcf_candidate_bases_idx.sort(key=lambda x: x != idx_of_ref_base)

            # Convert the genotypes to standard format
            # and count the number of each non-missing genotype
            genotype_counts = []
            for idx, idxval in enumerate(vcf_candidate_bases_idx):
                print(idx, idxval)
                fields[9:] = [f"{idx}/{idx}" if x == str(idxval) else x for x in fields[9:]]
                genotype_counts.append(fields[9:].count(f'{idx}/{idx}'))

            # New list of vcf bases after re-ordered and reference base identified
            vcf_final_bases = [vcf_candidate_bases[i] for i in vcf_candidate_bases_idx]

            # Replace the old vcf reference and alternative bases with the new ones
            fields[3] = vcf_final_bases[0]
            vcf_final_bases.pop(0)
            fields[4] = ",".join(vcf_final_bases)

            # Save reference counts and frequencies
            total_genotype_count = sum(genotype_counts)
            ref_counts_value = genotype_counts.pop(0)
            ref_counts = f"RC={ref_counts_value}"
            ref_frequencies = f"RF={round(ref_counts_value/total_genotype_count, ndigits=4)}"
            
            # Save alternative counts and frequencies
            alt_counts = f"AC={','.join(str(x) for x in genotype_counts)}"
            alt_frequencies = f"AF={','.join(str(round(x / sum(genotype_counts), ndigits=4)) for x in genotype_counts)}"
            fields[7] = alt_counts + ";" + alt_frequencies

            # Write the new line
            reassembled_line = '\t'.join(fields)

            # Write the new line
            with open(output_file, 'a', encoding='utf-8') as output_file_obj:
                output_file_obj.write(reassembled_line + "\n")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser("python remake_vcf.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input vcf file", dest="input_file")
    parser.add_argument("-r", "--reference_file", type=str, required=True, help="Reference file", dest="reference_file")
    parser.add_argument("-c", "--chrom_name", type=str, required=True, help="Chromosome name", dest="chrom_name")
    parser.add_argument("-l", "--chrom_length", type=int, required=True, help="Chromosome length", dest="chrom_length")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Output vcf file", dest="output_file")
    parser.add_argument("-s", "--samtools_path", type=str, required=True, help="Path to samtools", dest="samtools_path")
    return parser


def main(argv: List[str]) -> None:
    """
    Main function
    """
    # Parse command line arguments
    parser = parse_args()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define input and output files
    input_file = Path(args.input_file)
    reference_file = Path(args.reference_file)
    chrom_name = args.chrom_name
    chrom_length = args.chrom_length
    output_file = Path(args.output_file)
    samtools_path = Path(args.samtools_path)

    # # Define input and output files
    # input_file = Path("remake_vcf/example/example.vcf")
    # reference_file = Path("remake_vcf/example/reference/dm3.fa")
    # chrom_name = "chr2L"
    # chrom_length = 23011544
    # output_file = Path("remake_vcf/example/example_output.vcf")
    # samtools_path = Path("/usr/local/anaconda3/envs/bioinfo/bin/samtools")

    # Run remake_vcf
    remake_vcf(input_file, reference_file, chrom_name, chrom_length, output_file, samtools_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
