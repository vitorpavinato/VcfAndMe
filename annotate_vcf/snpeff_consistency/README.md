# SNPEff Consistency

A tool for analyzing multiple annotations from SNPEff output, providing both strict and majority-rule based approaches for effect classification.

## Overview

SNPEff often assigns multiple annotations to a single SNP, especially in compact genomes with overlapping genomic elements. This tool provides flexible approaches to handle such multiple annotations:

1. **Strict Effect Analysis**: Ideal for sparse genomes (e.g., human) where SNPs often have single, clear effects
2. **Majority Rule Analysis**: Better suited for compact genomes (e.g., Drosophila) with frequent overlapping features
3. **Custom Annotation Integration**: Currently supports short-intron annotations, with plans for expansion

## Usage

```zsh
python snpeff_consistency.py input.vcf -o output.txt -t 0.75
```

## Arguments

- `input.vcf`: Input VCF file with SNPEff annotations
- `-o, --output`: Output file name (default: annotation_summary.txt)
- `-t, --threshold`: Threshold for majority rule (default: 0.75)

## Output Format
The tool produces a tab-delimited file with the following columns:

- `chrom`: Chromosome
- `pos`: Position
- `strict_effect_bool`: True if SNP has only one effect
- `strict_effect_name`: The effect name (or "undefined" if multiple effects)
- `majority_rule_bool`: True if any effect appears above threshold
- `majority_rule_effect`: The majority effect, with special cases:
    - `Effect name`: When a clear majority exists
    - `tie:EFFECT1=EFFECT2`: When two effects have equal counts
    - `unclear`:EFFECT=XX%: When no effect reaches threshold
    - `+SI` suffix: Indicates presence of short intron annotation

## Features
- Handles position-based effects (Upstream/Downstream/UTRs) with distance information
- Integrates custom annotations (currently short-introns)
- Provides clear indication of ambiguous cases
- Configurable majority threshold

## Future Development
### Planned Features

- Distance-based majority rule implementation
- Abstract handling of custom annotations
- Support for multiple custom annotation types
- Enhanced distance-based analysis

## Implementation Details
The distance-based majority rule will consider both effect type and distance to features, allowing for more nuanced analysis of position-dependent effects.

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Citation
[Add citation information if applicable]