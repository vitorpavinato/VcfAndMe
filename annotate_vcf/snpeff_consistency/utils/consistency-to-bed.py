"""
Convert the consistency table to a BED file
"""

#!/usr/bin/env python3

import argparse

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Convert filtered consistency table to BED format')
    
    parser.add_argument('input',
                       help='Input filtered table from filter_consistency.py')
    
    parser.add_argument('-o', '--output',
                       default='consistency.bed',
                       help='Output BED file name (default: consistency.bed)')
    
    return parser.parse_args()

def table_to_bed(input_file, output_file):
    """Convert filtered consistency table to BED format"""
    with open(input_file, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:

        # Write header required by vcftools
        outfile.write("chr\tstart\tend\teffect\textra_info\n")

        # Skip header
        _ = infile.readline()

        # Store entries to sort
        entries = []

        for line in infile:
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])

            # Get effect type and other fields based on mode
            # (we'll take whatever effect column is present)
            effect = None
            for field in fields[2:]:
                if field in ['INTRON', 'SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING']:
                    effect = field
                    break

            if effect is None:
                continue

            # BED format is 0-based, half-open
            bed_start = pos - 1
            bed_end = pos

            # Include all original fields as extra columns
            extra_fields = '::'.join(fields[2:])

            entries.append((
                chrom,
                bed_start,
                bed_end,
                effect,
                extra_fields
            ))

        # Sort entries by chromosome and position
        entries.sort(key=lambda x: (x[0], x[1]))

        # Write sorted entries
        for entry in entries:
            outfile.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\n")

def main():
    """Main function"""
    args = parse_arguments()
    table_to_bed(args.input, args.output)

if __name__ == "__main__":
    main()
