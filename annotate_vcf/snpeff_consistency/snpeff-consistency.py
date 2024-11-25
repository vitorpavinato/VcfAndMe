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

import argparse
from collections import defaultdict, Counter
import re
from typing import Dict, List, Set, Tuple


class EffectAnalyzer:
    """
    Analyzer for SNPEff annotations with support for standard and custom effects.
    
    Parameters:
        majority_threshold (float): Minimum proportion to consider an effect as majority (default: 0.75)
    """
    
    def __init__(self, majority_threshold: float = 0.75):
        self.majority_threshold = majority_threshold
        self.results = []
        self.reset_statistics()
        
    def reset_statistics(self):
        """Reset all statistical counters"""
        self.total_snps = 0
        self.conflicting_snps = 0
        self.strict_consistent_snps = 0
        self.snps_with_short_introns = 0
        self.stats = defaultdict(int)

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
            return {'type': effect_type}
            
        details = detail_match.group(1).split('|')
        
        # Base parsed information
        parsed = {
            'type': effect_type,
            'impact': details[0] if len(details) > 0 else None,
            'transcript': details[8] if len(details) > 8 else None,
            'gene': details[5] if len(details) > 5 else None,
            'feature_type': details[6] if len(details) > 6 else None,
            'coding_status': details[7] if len(details) > 7 else None
        }
        
        # Add effect-specific details
        if effect_type in ['UPSTREAM', 'DOWNSTREAM', 'UTR_5_PRIME', 'UTR_3_PRIME']:
            parsed.update({
                'distance': int(details[2]) if details[2].isdigit() else None,
                'detail_type': 'position'
            })
            
        elif effect_type in ['SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING', 'STOP_GAINED', 'STOP_LOST', 'START_LOST', 'SYNONYMOUS_STOP', 'NON_SYNONYMOUS_START']:
            parsed.update({
                'codon_change': details[2] if len(details) > 2 else None,
                'aa_position': details[3] if len(details) > 3 else None,
                'cds_size': int(details[4]) if len(details) > 4 and details[4].isdigit() else None,
                'exon_rank': int(details[9]) if len(details) > 9 and details[9].isdigit() else None,
                'detail_type': 'coding'
            })
            
        elif effect_type == 'INTRON':
            parsed.update({
                'intron_rank': int(details[9]) if len(details) > 9 and details[9].isdigit() else None,
                'detail_type': 'intron'
            })

        elif effect_type == 'EXON':
            parsed.update({
                'exon_rank': int(details[9]) if len(details) > 9 and details[9].isdigit() else None,
                'detail_type': 'exon'
            })
            
        return parsed

    def _parse_custom_effect(self, effect_string: str) -> Dict:
        """
        Parse custom annotation effects (e.g., short introns).
        """
        custom_type = re.search(r'CUSTOM\[(.*?)\]', effect_string).group(1)
        
        detail_match = re.search(r'\((.*?)\)', effect_string)
        if not detail_match:
            return {'type': 'CUSTOM', 'custom_type': custom_type}
            
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
        effects = [self.parse_effect_detail(e) for e in effects_raw]
        
        # Group effects by their detail type
        effects_by_type = defaultdict(list)
        for effect in effects:
            effects_by_type[effect.get('detail_type', 'other')].append(effect)

        # Analyze custom annotations
        custom_annotations = self._analyze_custom_annotations(effects_by_type.get('custom', []))
        
        result = {
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt,
            'effects_by_type': dict(effects_by_type),
            'custom_annotations': custom_annotations
        }
        
        self.results.append(result)
        return result

    def _analyze_custom_annotations(self, custom_effects: List[Dict]) -> Dict:
        """
        Analyze custom annotations (e.g., short introns).
        """
        if not custom_effects:
            return {'has_custom_annotations': False}
            
        custom_by_type = defaultdict(list)
        for effect in custom_effects:
            custom_by_type[effect['custom_type']].append(effect)
        
        # Analyze short introns if present
        short_intron_analysis = None
        if 'dm6_short_introns' in custom_by_type:
            short_intron_effects = custom_by_type['dm6_short_introns']
            short_intron_analysis = {
                'total_annotations': len(short_intron_effects),
                'unique_transcripts': len({e['transcript'] for e in short_intron_effects if e['transcript']}),
                'intron_numbers': sorted(list({e['intron_number'] for e in short_intron_effects if e['intron_number'] is not None}))
            }
            
            self.snps_with_short_introns += 1
        
        return {
            'has_custom_annotations': True,
            'annotation_types': list(custom_by_type.keys()),
            'short_intron_analysis': short_intron_analysis
        }

    def create_simple_summary_table(self) -> str:
        """
        Create summary table with separate columns for boolean consistency checks and effect names.
        """
        output = "chrom\tpos\tstrict_effect_bool\tstrict_effect_name\tmajority_rule_bool\tmajority_rule_effect\n"
        
        for result in self.results:
            chrom = result['chrom']
            pos = result['pos']
            
            # Separate standard and custom effects
            standard_effects = []
            for effect_type, effects in result['effects_by_type'].items():
                if effect_type != 'custom':  # Only consider standard annotations
                    standard_effects.extend(effects)
            
            # Count standard effects
            effect_counts = Counter(e['type'] for e in standard_effects)
            total_effects = sum(effect_counts.values())
            
            # Strict effect analysis
            strict_consistent = len(effect_counts) == 1
            strict_effect_name = list(effect_counts.keys())[0] if strict_consistent else "undefined"
            
            # Majority rule analysis
            majority_effect = "undefined"
            has_majority = False
            
            if effect_counts:
                # Get top two most common effects
                top_effects = effect_counts.most_common(2)
                most_common_effect = top_effects[0]
                effect_proportion = most_common_effect[1] / total_effects
                
                # Check if we have a clear majority
                if effect_proportion >= self.majority_threshold:
                    has_majority = True
                    majority_effect = most_common_effect[0]
                elif len(top_effects) > 1:
                    # Check if we have a tie or very close proportions
                    second_effect = top_effects[1]
                    if most_common_effect[1] == second_effect[1]:
                        majority_effect = f"tie:{most_common_effect[0]}={second_effect[0]}"
                    else:
                        majority_effect = f"unclear:{most_common_effect[0]}={effect_proportion:.2%}"
            
            # Handle custom annotations (short introns)
            if result['custom_annotations']['has_custom_annotations']:
                if result['custom_annotations']['short_intron_analysis']:
                    if strict_effect_name != "undefined":
                        strict_effect_name += "+SI"
                    if majority_effect != "undefined" and not majority_effect.startswith(("tie:", "unclear:")):
                        majority_effect += "+SI"
            
            # Build output line
            output += f"{chrom}\t{pos}\t{str(strict_consistent)}\t{strict_effect_name}\t{str(has_majority)}\t{majority_effect}\n"
        
        return output


def parse_arguments():
   """Parse command line arguments"""
   
   # Parse command line arguments
   parser = argparse.ArgumentParser(description='Analyze SNPEff annotations consistency in VCF file')
   
   parser.add_argument('-i', '--input', help='Input VCF file')
   
   parser.add_argument('-o', '--output', 
                      default='annotation_summary.txt',
                      help='Output file name (default: annotation_summary.txt)')
   
   parser.add_argument('-t', '--threshold',
                      type=float,
                      default=0.75,
                      help='Threshold for majority rule (default: 0.75)')
   
   return parser.parse_args()


def main():
   """Main function"""
   args = parse_arguments()
   
   # Initialize analyzer
   analyzer = EffectAnalyzer(majority_threshold=args.threshold)
   
   # Process VCF file
   with open(args.input, 'r') as vcf:
       for line in vcf:
           if line and not line.startswith('#'):
               analyzer.analyze_effects(line)
   
   # Create and save summary table
   summary_table = analyzer.create_simple_summary_table()
   with open(args.output, 'w') as f:
       f.write(summary_table)


if __name__ == "__main__":
    main()


