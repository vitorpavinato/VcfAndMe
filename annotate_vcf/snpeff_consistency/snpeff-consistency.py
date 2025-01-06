"""
This module provides a class for analyzing SNPEff annotations in VCF files
particularly to determine if the annotations are consistent with each other.

Script renamed from vcf_annotation_analyzer_v4.py (or snpeff-annotation-analyzer.py from Claude.ai)

SNPEff Annotation Analyzer
-------------------------
A module for analyzing SNPEff annotations with support for:
- Effect-specific parsing
- Standard and custom annotation handling
- Strict and majority rule consistency analysis
"""

import os
import argparse
import logging
from collections import Counter
import re
from typing import Dict, List, Tuple


class EffectAnalyzer:
    """
    Analyzer for SNPEff annotations with support for standard and custom effects.
    
    Parameters:
        majority_threshold (float): Minimum proportion to consider an effect as majority (default: None)
        distance_threshold (fload): Minimum distance between a SNP and a feature to evaluate the effect (default: None)
    """
 
    # Define effect types - hard-coded
    # Positional effects
    POSITION_BASED_EFFECTS = {
        'DOWNSTREAM',
        'UPSTREAM', 
        'UTR_3_PRIME', 
        'UTR_5_PRIME', 
        'EXON', #  moved from feature-based
        'INTERGENIC', #  moved from feature-based
        'INTRON', #  moved from feature-based
        'SPLICE_SITE_REGION', #  moved from feature-based
        'SPLICE_SITE_REGION+EXON', #  moved from feature-based
        'SPLICE_SITE_REGION+INTRON' #  moved from feature-based
    }

    # Feature-based effects
    FEATURE_BASED_EFFECTS = {
        'NON_SYNONYMOUS_CODING',
        'NON_SYNONYMOUS_CODING+SPLICE_SITE_REGION',
        'NON_SYNONYMOUS_START',
        'SPLICE_SITE_ACCEPTOR+INTRON',
        'SPLICE_SITE_DONOR+INTRON',
        'SPLICE_SITE_REGION+SYNONYMOUS_CODING',
        'SPLICE_SITE_REGION+SYNONYMOUS_STOP',
        'START_GAINED',
        'START_LOST',
        'STOP_GAINED',
        'STOP_LOST',
        'STOP_LOST+SPLICE_SITE_REGION',
        'SYNONYMOUS_CODING',
        'SYNONYMOUS_STOP'
    }

    SPECIFIC_EFFECTS = {
        'INTERGENIC',
        'INTRON',
        'SYNONYMOUS_CODING',
        'NON_SYNONYMOUS_CODING'
    }
    
    def __init__(self, mode: str, majority_threshold: float = None, distance_threshold: int = None):
        self.mode = mode
        self.majority_threshold = majority_threshold
        self.distance_threshold = distance_threshold
        self.results = []
        self.reset_statistics()
        
    def reset_statistics(self):
        """Reset all statistical counters"""
        self.total_snps = 0
        self.consistent_snps = 0
        self.inconsistent_snps = 0
        self.distance_filtered_snps = 0
        self.snps_with_custom_annotations = 0

    def parse_effect_detail(self, effect_string: str) -> Dict:
        """
        Parse detailed information from an effect annotation.
        
        Args:
            effect_string: Raw effect string from VCF
            
        Returns:
            Dictionary containing parsed effect information
        """
        # Check if this is a custom annotation
        if effect_string.startswith('CUSTOM'):
            return self._parse_custom_effect(effect_string)
        
        effect_type = re.match(r'([^(]+)', effect_string).group(1)
        detail_match = re.search(r'\((.*?)\)', effect_string)
        
        if not detail_match:
            return {'type': effect_type, 'detail_type': 'standard'}
            
        details = detail_match.group(1).split('|')
        
        # Base parsed information
        parsed = {
            'type': effect_type,
            'impact': details[0] if len(details) > 0 else None,
            'transcript': details[8] if len(details) > 8 else None,
            'gene': details[5] if len(details) > 5 else None,
            'feature_type': details[6] if len(details) > 6 else None,
            'coding_status': details[7] if len(details) > 7 else None,
            'detail_type': 'standard'
        }
        
        # Add effect-specific details
        if effect_type in self.POSITION_BASED_EFFECTS:
            parsed.update({
                'distance': int(details[2]) if details[2].isdigit() else None,
                'detail_type': 'position'
            })
            
        elif effect_type in self.FEATURE_BASED_EFFECTS:
            parsed.update({
                'functional_class': details[1] if len(details) > 1 else None,
                'codon_change': details[2] if len(details) > 2 else None,
                'aa_position': details[3] if len(details) > 3 else None,
                'cds_size': int(details[4]) if len(details) > 4 and details[4].isdigit() else None,
                'feature_rank': int(details[9]) if len(details) > 9 and details[9].isdigit() else None,
                'detail_type': 'feature'
            })

        return parsed

    def _parse_custom_effect(self, effect_string: str) -> Dict:
        """
        Parse custom annotation effects (e.g., short introns).
        """
        custom_type = re.search(r'CUSTOM\[(.*?)\]', effect_string).group(1)
        
        detail_match = re.search(r'\((.*?)\)', effect_string)
        if not detail_match:
            return {
                'type': 'CUSTOM', 
                'custom_type': custom_type,
                'detail_type': 'custom'
            }
            
        details = detail_match.group(1).split('|')
        
        transcript_info = details[6] if len(details) > 5 else ''
        transcript_match = re.match(r'(NM_\d+\.\d+)_intron_(\d+)', transcript_info)
        
        return {
            'type': 'CUSTOM',
            'custom_type': custom_type,
            'impact': details[0] if len(details) > 0 else None,
            'transcript': transcript_match.group(1) if transcript_match else None,
            'intron_number': int(transcript_match.group(2)) if transcript_match else None,
            'raw_transcript_info': transcript_info,
            'detail_type': 'custom'
        }

    def analyze_effects(self, vcf_line: str) -> Dict:
        """
        Analyze effects from a VCF line.
        """
        # Skip header lines
        if vcf_line.startswith('#'):
            return None
            
        # Parse VCF fields
        fields = vcf_line.strip().split('\t')
        chrom, pos, _, ref, alt = fields[0:5]
        
        # Extract EFF field
        info = fields[7]
        eff_match = re.search(r'EFF=([^;]+)', info)
        if not eff_match:
            return None
            
        # Parse effects
        effects_raw = eff_match.group(1).split(',')
        all_effects = [self.parse_effect_detail(e) for e in effects_raw]

        # Separate standard and custom effects
        standard_effects = [e for e in all_effects if e['detail_type'] != 'custom']
        custom_effects = [e for e in all_effects if e['detail_type'] == 'custom']
        
        # Update statistics
        self.total_snps += 1
        if custom_effects:
            self.snps_with_custom_annotations += 1

        # Analyze based on mode
        if self.mode == 'strict':
            result = self._analyze_strict(standard_effects)
        elif self.mode == 'rule':
            if standard_effects and standard_effects[0]['type'] in self.POSITION_BASED_EFFECTS:
                result = self._analyze_majority_rule(standard_effects)
            else:
                result = self._analyze_first_effect(standard_effects)
        else:  # mode == 'specific'
            result = self._analyze_specific_effects(standard_effects)

        # Add custom annotation analysis
        custom_analysis = self._analyze_custom_effects(custom_effects)

        # Combine results
        combined_result = {
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt,
            'is_consistent': result['is_consistent'],
            'effect_type': result['effect_type'],
            'analysis_method': result['analysis_method'],
            'custom_annotations': custom_analysis,
            'all_effects': [e['type'] for e in standard_effects],
            'filtered_effects': result.get('filtered_effects', None)
        }
        
        # Update statistics
        if combined_result['is_consistent']:
            self.consistent_snps += 1
        else:
            self.inconsistent_snps += 1
        
        self.results.append(combined_result)
        return combined_result

    
    def _analyze_strict(self, effects: List[Dict]) -> Dict:
        """Analyze using strict consistency rule"""
        unique_effects = set(e['type'] for e in effects)
        return {
            'effect_type': list(unique_effects)[0] if len(unique_effects) == 1 else 'undefined',
            'is_consistent': len(unique_effects) == 1,
            'analysis_method': 'strict'
        }


    def _analyze_majority_rule(self, effects: List[Dict]) -> Dict:
        """Analyze using majority rule, properly handling both position and feature-based effects"""
        if self.distance_threshold is not None:
            filtered_effects = [
                effect for effect in effects
                if(effect['detail_type'] == 'feature' or #  Keep all feature-based effects
                   effect.get('distance') is None or #  Keep effects without distance information
                   effect['distance'] <= self.distance_threshold) #  Filter by distance threshold
            ]
            if len(filtered_effects) != len(effects):
                self.distance_filtered_snps += 1
        else:
            filtered_effects = effects
        
        effect_counts = Counter(e['type'] for e in filtered_effects)
        total_effects = len(filtered_effects)
        
        if total_effects == 0:
            return {
                'effect_type': 'undefined',
                'is_consistent': False,
                'analysis_method': 'majority_rule',
                'filtered_effects': []
            }
        
        majority_effect = effect_counts.most_common(1)[0]
        has_majority = majority_effect[1] / total_effects >= self.majority_threshold
        
        return {
            'effect_type': majority_effect[0] if has_majority else 'undefined',
            'is_consistent': has_majority,
            'analysis_method': 'majority_rule',
            'filtered_effects': [e['type'] for e in filtered_effects]
        }

    
    def _analyze_first_effect(self, effects: List[Dict]) -> Dict:
        """Analyze using first effect rule"""
        if not effects:
            return {
                'effect_type': 'undefined',
                'is_consistent': False,
                'analysis_method': 'first_effect'
            }
        
        return {
            'effect_type': effects[0]['type'],
            'is_consistent': True,
            'analysis_method': 'first_effect'
        }


    def _analyze_specific_effects(self, effects: List[Dict]) -> Dict:
        """Analyze looking for specific effects of interest"""
        found_effects = {e['type'] for e in effects if e['type'] in self.SPECIFIC_EFFECTS}
        
        if not found_effects:
            return {
                'effect_type': 'undefined',
                'is_consistent': False,
                'analysis_method': 'specific'
            }
        
        if len(found_effects) == 1:
            return {
                'effect_type': list(found_effects)[0],
                'is_consistent': True,
                'analysis_method': 'specific'
            }
        
        # Handle conflicts based on severity
        if 'NON_SYNONYMOUS_CODING' in found_effects:
            return {
                'effect_type': 'NON_SYNONYMOUS_CODING',
                'is_consistent': True,
                'analysis_method': 'specific',
                'note': 'chosen effect arbitrarily'
            }
        elif 'SYNONYMOUS_CODING' in found_effects:
            return {
                'effect_type': 'SYNONYMOUS_CODING',
                'is_consistent': True,
                'analysis_method': 'specific',
                'note': 'chosen effect arbitrarily'
            }
        elif 'INTRON' in found_effects:
            return {
                'effect_type': 'INTRON',
                'is_consistent': True,
                'analysis_method': 'specific',
                'note': 'chosen effect arbitrarily'
            }
        else:
            return {
                'effect_type': 'INTERGENIC',
                'is_consistent': True,
                'analysis_method': 'specific',
                'note': 'chosen effect arbitrarily'
            }

    def _analyze_custom_effects(self, effects: List[Dict]) -> Dict:
        """Analyze custom annotations"""
        if not effects:
            return {
                'has_custom_annotations': False,
                'annotation_types': []
            }
        
        custom_types = {e['custom_type'] for e in effects if e.get('custom_type')}
        
        return {
            'has_custom_annotations': True,
            'annotation_types': list(custom_types)
        }

    
    def create_summary_table(self) -> str:
        """Create summary table with separate custom annotation columns"""
        
        # Base headers for all modes
        base_header = "chrom\tpos\t"
        custom_header = "has_custom_annotation\tcustom_annotation_type\t"

        # Mode-specific headers
        mode_header = None
        if self.mode == 'strict':
            mode_header = "strict_effect_bool\tstrict_effect_name"
        elif self.mode == 'rule':
            mode_header = "rule_effect_bool\trule_effect_name\tdistance_filtered"
        elif self.mode == 'specific':
            mode_header = "specific_effect_bool\tspecific_effect_name"

        output = base_header + custom_header + mode_header + "\n"

        for result in self.results:
            # Base info
            chrom = result['chrom']
            pos = result['pos']
            base_info = f"{chrom}\t{pos}\t"
            
            # Custom annotation info
            has_custom = result['custom_annotations']['has_custom_annotations']
            # custom_type = result['custom_annotations']['annotation_types'][0] if has_custom else 'NA'
            custom_types = '+'.join(result['custom_annotations']['annotation_types']) if has_custom else 'NA'
            custom_info = f"{str(has_custom)}\t{custom_types}\t"
            
            # Mode-specific info
            is_consistent = result['is_consistent']
            effect_type = result['effect_type']

            if self.mode == 'strict':
                mode_info = f"{str(is_consistent)}\t{effect_type}"

            elif self.mode == 'rule':
                distance_filtered = "NA"  # this initializes the variable
                if result['analysis_method'] == 'majority_rule':  # Only for position-based effects
                    filtered_effects = result.get('filtered_effects', [])  # Default to empty
                    all_effects = result.get('all_effects', [])
                    
                    # Check if any effects were filtered out
                    if len(filtered_effects) != len(all_effects):
                        distance_filtered = "True"
                    else:
                        distance_filtered = "False"

                else:  # first-effect method
                    distance_filtered = "not_applicable"

                mode_info = f"{str(is_consistent)}\t{effect_type}\t{distance_filtered}"

            else:  # mode == 'specific'
                # conflict_resolution = result.get('note', 'NA')
                mode_info = f"{str(is_consistent)}\t{effect_type}"

            output += base_info + custom_info + mode_info + "\n"

        return output


    def get_summary_stats(self) -> Dict:
        """Get analysis summary statistics"""

        summary = {
            'total_snps': self.total_snps,
            'consistent_snps': self.consistent_snps,
            'inconsistent_snps': self.inconsistent_snps,
            'snps_with_custom_annotations': self.snps_with_custom_annotations
        }

        if self.mode == 'rule':
            summary['distance_filtered_snps'] = self.distance_filtered_snps

        # Update the log with the summary of the analysis
        logging.info("Analysis Completed Successfully!")
        logging.info("Summary of the Analysis:")
        logging.info("Total SNPs processed: %s", summary['total_snps'])
        logging.info("Consistent annotations: %s", summary['consistent_snps'])
        logging.info("Inconsistent annotations: %s", summary['inconsistent_snps'])
        logging.info("SNPs with custom annotations: %s", summary['snps_with_custom_annotations'])

        if self.mode == 'rule':
            logging.info("SNPs with distance filtering: %s", summary['distance_filtered_snps'])

        return summary

    def get_detailed_stats(self) -> Dict:
        """Get detailed counts of different annotations and effects"""
        
        summary = {
            'custom_annotations': Counter(
                result['custom_annotations']['has_custom_annotations']
                for result in self.results
            ),
            'custom_types': Counter(
                result['custom_annotations']['annotation_types'][0]
                if result['custom_annotations']['annotation_types']
                else 'NA'
                for result in self.results
            ),
            'effect_consistency': Counter(
                result['is_consistent']
                for result in self.results
            ),
            'effect_names': Counter(
                result['effect_type']
                for result in self.results
            )
        }

        # # Log the detailed counts
        # logging.info("Detailed Summary:")
        # logging.info("Custom Annotations:")
        # logging.info("Has custom annotations: %s", dict(summary['custom_annotations']))
        # logging.info("Custom annotation types: %s", dict(summary['custom_types']))
        # logging.info("Effect Summary:")
        # logging.info("Effect consistency: %s", dict(summary['effect_consistency']))
        # logging.info("Effect names distribution:")
        # for effect, count in summary['effect_names'].items():
        #     logging.info("%s: %s", effect, count)

        return summary


def validate_parameters(args):
    """Validate parameter combinations and ranges"""

    logging.basicConfig(
        filename='snpeff_consistency.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Check mode
    if args.mode not in ['strict', 'rule', 'specific']:
        raise ValueError(f"Invalid mode: {args.mode}")

    # Strict or specific mode
    if args.mode in ('strict', 'specific'):
        logging.info("Mode: %s", args.mode)

        if args.threshold != 0.75 or args.distance is not None:
            logging.warning("Threshold (-t) and distance (-d) parameters are ignored in strict mode")

    else: # Rule mode
        logging.info("Mode: %s", args.mode)

        if args.threshold is not None:
            if not 0 < args.threshold <= 1:
                raise ValueError(f"Threshold must be between 0 and 1, got {args.threshold}")

        if args.distance is not None:
            if args.distance <= 0:
                raise ValueError(f"Distance threshold must be positive, got {args.distance}")

        # Rule mode parameters for the majority rule
        logging.info(" + Majority rule threshold: %s", args.threshold)
        logging.info(" + Distance threshold: %s", args.distance if args.distance is not None else 'None (using all effects)')

    # Check if output file already exists
    if os.path.exists(f"{args.output}.txt") or os.path.exists(f"{args.output}_stats.txt"):
        logging.warning("Output files %s* already exist and will be overwritten", args.output)
        print(f"Warning: Output files {args.output}* already exist and will be overwritten")

    return args


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analyze SNPEff annotations in VCF file')

    parser.add_argument('input',
                       help='Input VCF file')

    parser.add_argument('-o', '--output',
                       default='annotation_summary',
                       help='Output file name (default: annotation_summary.txt)')

    parser.add_argument('--mode',
                       choices=['strict', 'rule', 'specific'],
                       required=True,
                       help='Analysis mode (strict: all effects must be the same; '
                            'rule: apply majority rule for position-based effects and '
                            'first effect for feature-based effects; '
                            'specific: look for specific effects of interest)')

    parser.add_argument('-t', '--threshold',
                       type=float,
                       default=0.75,
                       help='Threshold for majority rule when using position-based effects '
                            '(default: 0.75, only used in rule mode)')

    parser.add_argument('-d', '--distance',
                       type=int,
                       default=None,
                       help='Distance threshold in bp for position-based effects. '
                            'When set, only effects within this distance are considered '
                            'for majority rule (optional, only used in rule mode)')

    parser.add_argument('--stats',
                       action='store_true',
                       help='Print summary statistics (default: False)')

    args = parser.parse_args()
    return validate_parameters(args)

def main():
    """Main function"""
    args = parse_arguments()

    # Initialize analyzer
    analyzer = EffectAnalyzer(
        mode=args.mode,
        majority_threshold=args.threshold if args.mode == 'rule' else None,
        distance_threshold=args.distance if args.mode == 'rule' else None
    )

    # Process VCF file
    with open(args.input, 'r', encoding='utf-8') as vcf:
        for line in vcf:
            if line and not line.startswith('#'):
                analyzer.analyze_effects(line)

    # Create and save summary table
    summary_table = analyzer.create_summary_table()
    with open(f'{args.output}.txt', 'w', encoding='utf-8') as f:
        f.write(summary_table)

    # Get the simple summary
    summary_stats = analyzer.get_summary_stats()
    with open(f'{args.output}_stats.txt', 'w', encoding='utf-8') as f:
        f.write("=== Summary of the Analysis ===\n")
        f.write(f"Total SNPs processed: {summary_stats['total_snps']}\n")
        f.write(f"Consistent annotations: {summary_stats['consistent_snps']}\n")
        f.write(f"Inconsistent annotations: {summary_stats['inconsistent_snps']}\n")
        f.write(f"SNPs with custom annotations: {summary_stats['snps_with_custom_annotations']}\n")
        if args.mode == 'rule':
            f.write(f"SNPs with distance filtering: {summary_stats['distance_filtered_snps']}\n")
        f.write("\n")

    if args.stats:
        # Get detailed summary stats
        detailed_summary_stats = analyzer.get_detailed_stats()
        with open(f'{args.output}_stats.txt', 'a', encoding='utf-8') as f:
            f.write("=== Custom Annotation Stats ===\n")
            f.write(f"Has custom annotations: {dict(detailed_summary_stats['custom_annotations'])}\n")
            f.write(f"Custom annotation types: {dict(detailed_summary_stats['custom_types'])}\n\n")
            f.write("=== Effect Stats ===\n")
            f.write(f"Effect consistency: {dict(detailed_summary_stats['effect_consistency'])}\n\n")
            f.write("===Effect names distribution===\n")
            for effect, count in detailed_summary_stats['effect_names'].items():
                f.write(f"{effect}: {count}\n")


if __name__ == "__main__":
    main()
