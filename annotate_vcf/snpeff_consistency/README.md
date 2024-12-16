# SNPEff Consistency

A tool for analyzing and consolidating multiple SNPEff annotations in VCF files.

## Overview

SNPEff often assigns multiple annotations to a single SNP, especially in compact genomes with overlapping genomic elements.
SNPeff-consistency provides flexible approaches to handle multiple annotations that SNPEff assigns to a single SNP. It offers three analysis modes:
1. **Strict consistency checking**: Ideal for sparse genomes (e.g., human) where SNPs often have single, clear effects
2. **Rule-based analysis for position and feature-based effects**: Better suited for compact genomes (e.g., Drosophila) with frequent overlapping features
3. **Specific effect type identification**: This looks SNPs containing any of these effects: `NON_SYNONYMOUS_CODING`, `SYNONYMOUS_CODING`, `INTRONS`, and `INTERGENIC` (in this order)
4. **Custom Annotation Integration**: Currently supports short-intron annotations, with plans for expansion

## Features

- Multiple analysis modes for different use cases
- Distance-based filtering for positional effects (*e.g.* `DOWNSTREAM`, `UPSTREAM`, etc)
- Custom annotation handling (e.g., short introns)
- Detailed statistics and logging
- Conservative handling of functional annotations specially for specific type identification

## Installation

```bash
git clone [repository-url]
cd snpeff-consistency
```

## Dependencies
- Python 3.x
- No additional packages required

## Usage

Basic usage:
```zsh
python snpeff-consistency.py -o output --mode [strict|rule|specific] -t threshold -d distance --stats input.vcf
```

## Arguments

- `input.vcf`: Input VCF file with SNPEff annotations
- `--mode`: Analysis mode (required)
    - `strict`: All effects must be identical
    - `rule`: Apply majority rule for position-based effects and first effect for feature-based
    - `specific`: Look for specific effects of interest
- `-o, --output`: Output file name (default: annotation_summary.txt)
- `-t, --threshold`: Threshold for majority rule (default: 0.75)
- `-d, --distance`: Distance threshold in bp for position-based effects
- `--stats`:  Generate detailed statistics file

## Analysis Modes

### Strict Mode

Checks if all annotations for a SNP are identical. Useful for identifying SNPs with clear, unambiguous effects.

### Rule Mode

Applies different rules based on effect type:

- Position-based effects (UPSTREAM, DOWNSTREAM, etc.): Uses majority rule with optional distance filtering
- Feature-based effects (NON_SYNONYMOUS_CODING, etc.): Uses first effect (most severe)

### Specific Mode

Looks for specific effect types with priority handling:

- NON_SYNONYMOUS_CODING
- SYNONYMOUS_CODING
- INTERGENIC
- INTRON

## Output Files

1. Main output (tab-separeted):
```bash
chrom    pos    has_custom_annotation    custom_annotation_type    [mode-specific-columns]
```

2. Statistics file (when --stats is used):
- Custom annotation counts
- Effect type distribution
- Consistency statistics

## Exemples

1. Strict mode analysis:
```bash
python snpeff-consistency.py -o results --mode strict input.vcf 
```

2. Rule mode with distance threshold:
```bash
python snpeff-consistency.py -o results --mode rule -t 0.8 -d 500 input.vcf
```

3. Specific effect analysis with statistics:
```bash
python snpeff-consistency.py -o results.txt --mode specific --stats input.vcf
```

## Effect Types
### Position-based Effects

- DOWNSTREAM
- UPSTREAM
- UTR_3_PRIME
- UTR_5_PRIME
- EXON
- INTERGENIC
- INTRON
- SPLICE_SITE_REGION
- SPLICE_SITE_REGION+EXON
- SPLICE_SITE_REGION+INTRON

### Feature-based Effects

- NON_SYNONYMOUS_CODING
- SYNONYMOUS_CODING
- STOP_GAINED
- STOP_LOST
- START_LOST
And other coding/functional variants

## Logging
The tool creates a log file (snpeff_consistency.log) containing:

- Parameter settings
- Processing statistics
- Warnings and errors
- Analysis summary

## Future Development
### Planned Features

- Abstract handling of custom annotations
- Support for multiple custom annotation types

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Citation
[Add citation information if applicable]