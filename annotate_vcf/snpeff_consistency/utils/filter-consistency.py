"""
Filter SNPeff-consistency output based on mode and effects
"""

import argparse
import sys

# Define valid effects (from snpeff_consistency.py)
POSITION_BASED_EFFECTS = {
    'DOWNSTREAM',
    'UPSTREAM', 
    'UTR_3_PRIME', 
    'UTR_5_PRIME', 
    'EXON',
    'INTERGENIC',
    'INTRON',
    'SPLICE_SITE_REGION',
    'SPLICE_SITE_REGION+EXON',
    'SPLICE_SITE_REGION+INTRON'
}

FEATURE_BASED_EFFECTS = {
    'NON_SYNONYMOUS_CODING',
    'NON_SYNONYMOUS_START',
    'NON_SYNONYMOUS_CODING+SPLICE_SITE_REGION',
    'NON_SYNONYMOUS_START+SPLICE_SITE_REGION'
    'SPLICE_SITE_ACCEPTOR+INTRON',
    'SPLICE_SITE_DONOR+INTRON',
    'SPLICE_SITE_REGION+SYNONYMOUS_CODING',
    'SPLICE_SITE_REGION+SYNONYMOUS_STOP',
    'START_GAINED',
    'START_LOST',
    'START_LOST+SPLICE_SITE_REGION',
    'STOP_GAINED',
    'STOP_LOST',
    'STOP_LOST+SPLICE_SITE_REGION',
    'STOP_GAINED+SPLICE_SITE_REGION',
    'SYNONYMOUS_CODING',
    'SYNONYMOUS_STOP'
}

SPECIFIC_EFFECTS = {
        'INTERGENIC',
        'INTRON',
        'SYNONYMOUS_CODING',
        'NON_SYNONYMOUS_CODING'
    }

VALID_EFFECTS = set(POSITION_BASED_EFFECTS | FEATURE_BASED_EFFECTS | SPECIFIC_EFFECTS)

def validate_effects(effects):
    """
    Validate that all provided effects are valid
    """
    if not effects:
        return True
   
    invalid_effects = [effect for effect in effects if effect not in VALID_EFFECTS]
    if invalid_effects:
        print(f"Error: Invalid effects provided: {', '.join(invalid_effects)}", file=sys.stderr)
        print(f"Valid effects are: {', '.join(sorted(VALID_EFFECTS))}", file=sys.stderr)
        return False
    return True


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Filter SNPeff-consistency output based on mode and effects')

    parser.add_argument('input',
                       help='Input table from SNPeff-consistency')

    parser.add_argument('-o', '--output',
                       default='filtered_consistency.txt',
                       help='Output filtered table (default: filtered_consistency.txt)')

    parser.add_argument('--mode',
                       choices=['rule', 'specific', 'strict'],
                       required=True,
                       help='Analysis mode to filter')

    parser.add_argument('--effects',
                       nargs='+',
                       required=True,
                       help='List of effects to filter for')

    args = parser.parse_args()

    if not validate_effects(args.effects):
        sys.exit(1)

    return args

def filter_consistency(input_file, output_file, mode, effects):
    """
    Filter consistency table based on mode and effects
    """

    # Determine which columns to check based on mode
    if mode == 'rule':
        bool_col = 'rule_effect_bool'
        effect_col = 'rule_effect_name'
    elif mode == 'specific':
        bool_col = 'specific_effect_bool'
        effect_col = 'specific_effect_name'
    else:  # strict
        bool_col = 'strict_effect_bool'
        effect_col = 'strict_effect_name'

    with open(input_file, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
        # Get header
        header = infile.readline().strip()
        outfile.write(header + '\n')

        # Get column indices
        headers = header.split('\t')
        bool_idx = headers.index(bool_col)
        effect_idx = headers.index(effect_col)

        # Filter lines
        for line in infile:
            fields = line.strip().split('\t')
            # Check if effect matches and is consistent
            if (fields[effect_idx] in effects and fields[bool_idx].lower() == 'true'):
                outfile.write(line)

def main():
    """
    Main function
    """
    args = parse_arguments()
    filter_consistency(args.input, args.output, args.mode, args.effects)

if __name__ == "__main__":
    main()
