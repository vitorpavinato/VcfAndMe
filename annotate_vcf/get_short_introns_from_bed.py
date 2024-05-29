"""
This program will take a .BED file containing introns 
(downloaded from UCSC for example) it will filter out long 
introns retaining only introns shorter by the value defined in 
the short_intron_size parameter. This program also allows to 
trimm both end of the intron by the value defined in the 
trailling_size parameter.

This is a CLI of create_sfss/snp_utils.py filter_short_intron_from_bed().
Here it was renamed to get_short_introns_from_bed()
"""

import argparse
import sys


def get_short_introns_from_bed(
    inputfile: str, outputfile: str, chrom_list: list[str],
    short_intron_size: int = 86, trailling_size: int = 8
) -> None:
    """
    Function implementation that creates a .BED files with 
    only short introns intervals.
    """

    short_introns = []

    with open(inputfile, "r", encoding="utf-8") as inbed, open(outputfile, "w", encoding="utf-8") as outbed:
        for linen, line in enumerate(inbed):
            line = line.rstrip("\n")
            if (line.startswith("##") or line.startswith("#") or line.startswith(">") or
                    line.startswith("A") or line.startswith("T") or line.startswith("C") or line.startswith("G")):
                continue

            fields = [splits for splits in line.split("\t") if splits != ""]

            if len(fields) <= 1:
                continue

            # Unpack FIELDS items
            chrom, estart, eend, *_ = fields

            # Cast estarts and eends to integer
            estart = int(estart)
            eend = int(eend)

            # Process only chroms in the list:
            # This make the inclusion of any chromosome explicit.
            if any(x == chrom for x in chrom_list):
                if eend - estart < short_intron_size:
                    estart_trimmed = estart + trailling_size
                    eend_trimmed = eend - trailling_size
                    short_intron = [chrom, (estart_trimmed-1), eend_trimmed]
                    short_intron.extend(_)
                    outbed.write('\t'.join(str(item) for item in short_intron) + "\n")
                    short_introns.append(short_intron)

    return f"Short introns saved in: {outputfile}"


def parse_chromosome(chrom):
    """
    Parse chromsomes type
    """
    if chrom.startswith("chr"):
        return chrom
    else:
        return f"chr{chrom}"


def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser("python get_short_introns_from_bed.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="Input BED file name", dest="inputfile", required=True, type=str)
    parser.add_argument("-o", help="Output BED file name", dest="outputfile", required=True, type=str)
    parser.add_argument("-c", help="List of chroms to process", dest="chrom_list", nargs="+", required=True, type=parse_chromosome)
    parser.add_argument("-s", help="Short intron size", dest="short_intron_size", default=86, type=int)
    parser.add_argument("-t", help="Trailling size", dest="trailling_size", default=8, type=int)
    return parser


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
    chrom_list = args.chrom_list
    short_intron_size = args.short_intron_size
    trailling_size = args.trailling_size

    # Execute the function
    result = get_short_introns_from_bed(inputfile=inputfile,
                                        outputfile=outputfile,
                                        chrom_list=chrom_list,
                                        short_intron_size=short_intron_size,
                                        trailling_size=trailling_size)

    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
