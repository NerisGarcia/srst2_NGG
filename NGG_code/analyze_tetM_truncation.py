#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict
import re

# Add SRST2 script path to import the exact function
sys.path.append('/home/ngarcia/RESEARCH/Projects/2024_AMR_Enteroccoco/Goals/srst2_NGG/scripts/')

class TetMTruncationAnalyzer:
    def __init__(self, results_file, pileup_file):
        self.results_file = results_file
        self.pileup_file = pileup_file
        
        # SRST2 truncation detection parameters (exact from your srst2.py)
        self.HIGH_END_FRAC_THRESHOLD = 0.75   # 75% minimum threshold
        self.VERY_HIGH_END_FRAC = 0.90        # 90% - no coverage drop needed
        self.COVERAGE_DROP_RATIO = 3.0        # 3-fold coverage drop required
        self.MIN_DEPTH = 10                   # Minimum depth required
        self.EDGE_PAD = 10                    # Ignore positions near edges
        self.WINDOW_SIZE = 10                 # Look ahead window for coverage drop
        
        self.tetM_info = {}
        self.tetM_pileup_data = []
    
    def extract_tetM_from_results(self):
        """Extract tetM gene information from SRST2 results file"""
        print(f"Reading SRST2 results from: {self.results_file}")
        
        tetM_found = False
        
        with open(self.results_file, 'r') as f:
            header = f.readline().strip().split('\t')
            print(f"Header columns: {header}")
            
            # Find relevant column indices
            try:
                if 'gene' in header:  # fullgenes format
                    gene_idx = header.index('gene')
                    allele_idx = header.index('allele')
                    length_idx = header.index('length')
                    coverage_idx = header.index('coverage')
                    depth_idx = header.index('depth')
                    clusterid_idx = header.index('clusterid') if 'clusterid' in header else None
                    seqid_idx = header.index('seqid') if 'seqid' in header else None
                    
                    print(f"Looking for tetM in gene column (index {gene_idx})...")
                    
                    for line_num, line in enumerate(f, 2):
                        parts = line.strip().split('\t')
                        if len(parts) <= gene_idx:
                            continue
                            
                        gene = parts[gene_idx]
                        print(f"Line {line_num}: Found gene '{gene}'")
                        
                        # More flexible tetM detection
                        if any(tetm_variant in gene.lower() for tetm_variant in ['tetm', 'tet(m)', 'tet_m']):
                            tetM_found = True
                            self.tetM_info = {
                                'gene': gene,
                                'allele': parts[allele_idx] if allele_idx < len(parts) else 'N/A',
                                'length': int(parts[length_idx]) if length_idx < len(parts) else 0,
                                'coverage': float(parts[coverage_idx]) if coverage_idx < len(parts) else 0.0,
                                'depth': float(parts[depth_idx]) if depth_idx < len(parts) else 0.0,
                                'cluster_id': parts[clusterid_idx] if clusterid_idx and clusterid_idx < len(parts) else None,
                                'seq_id': parts[seqid_idx] if seqid_idx and seqid_idx < len(parts) else None
                            }
                            print(f"✅ Found tetM gene:")
                            print(f"  Gene: {self.tetM_info['gene']}")
                            print(f"  Allele: {self.tetM_info['allele']}")
                            print(f"  Length: {self.tetM_info['length']} bp")
                            print(f"  Coverage: {self.tetM_info['coverage']:.2f}%")
                            print(f"  Depth: {self.tetM_info['depth']:.2f}x")
                            if self.tetM_info['cluster_id']:
                                print(f"  Cluster ID: {self.tetM_info['cluster_id']}")
                            if self.tetM_info['seq_id']:
                                print(f"  Seq ID: {self.tetM_info['seq_id']}")
                            break
                            
            except ValueError as e:
                print(f"Error parsing results file: {e}")
                return False
        
        if not tetM_found:
            print("❌ No tetM gene found in SRST2 results")
            # Show all genes found for debugging
            print("Genes found in results file:")
            with open(self.results_file, 'r') as f:
                header = f.readline().strip().split('\t')
                gene_idx = header.index('gene') if 'gene' in header else 0
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > gene_idx:
                        print(f"  - {parts[gene_idx]}")
            return False
            
        return True
    
    def extract_tetM_pileup_data(self):
        """Extract pileup data for the tetM gene/allele"""
        print(f"\nExtracting tetM pileup data from: {self.pileup_file}")
        
        if not self.tetM_info:
            print("❌ No tetM info available - run extract_tetM_from_results first")
            return False
        
        # Build comprehensive list of possible reference names
        target_refs = []
        
        # Add all possible identifiers
        if self.tetM_info['cluster_id']:
            target_refs.append(self.tetM_info['cluster_id'])
        if self.tetM_info['seq_id']:
            target_refs.append(self.tetM_info['seq_id'])
            
        target_refs.extend([
            self.tetM_info['gene'],
            self.tetM_info['allele'],
            'tetM',
            'tetM_tetM',
            'tet(M)',
            'tetM_tetM-2'  # Based on your results showing tetM_tetM-2
        ])
        
        print(f"Searching for these reference names in pileup:")
        for target in target_refs:
            print(f"  - '{target}'")
        
        # First pass: find what references actually exist in pileup
        refs_in_pileup = set()
        with open(self.pileup_file, 'r') as f:
            for line in f:
                if line.strip():
                    ref_name = line.split('\t')[0]
                    refs_in_pileup.add(ref_name)
        
        print(f"\nReferences found in pileup file ({len(refs_in_pileup)}):")
        for ref in sorted(list(refs_in_pileup)[:20]):  # Show first 20
            print(f"  - {ref}")
        if len(refs_in_pileup) > 20:
            print(f"  ... and {len(refs_in_pileup) - 20} more")
        
        # Find matching references
        matching_refs = []
        for ref in refs_in_pileup:
            for target in target_refs:
                if target.lower() in ref.lower() or ref.lower() in target.lower():
                    if ref not in matching_refs:
                        matching_refs.append(ref)
                        print(f"✅ Match found: '{ref}' matches '{target}'")
        
        if not matching_refs:
            print("❌ No matching references found in pileup")
            # Try partial matches with 'tet'
            print("Trying partial matches with 'tet':")
            for ref in refs_in_pileup:
                if 'tet' in ref.lower():
                    matching_refs.append(ref)
                    print(f"  Partial match: {ref}")
        
        if not matching_refs:
            print("❌ Still no matches found")
            return False
        
        # Extract pileup data for matching references
        print(f"\nExtracting pileup data for {len(matching_refs)} matching reference(s)...")
        
        with open(self.pileup_file, 'r') as f:
            for line in f:
                if line.strip():
                    ref_name = line.split('\t')[0]
                    
                    if ref_name in matching_refs:
                        parts = line.strip().split('\t')
                        if len(parts) >= 5:
                            try:
                                pos = int(parts[1])
                                ref_base = parts[2]
                                depth = int(parts[3])
                                bases = parts[4] if len(parts) > 4 else ""
                                
                                self.tetM_pileup_data.append({
                                    'ref_name': ref_name,
                                    'pos': pos,
                                    'ref_base': ref_base,
                                    'depth': depth,
                                    'bases': bases
                                })
                            except (ValueError, IndexError):
                                continue
        
        if not self.tetM_pileup_data:
            print("❌ No pileup data extracted")
            return False
        
        # Sort by reference name and position
        self.tetM_pileup_data.sort(key=lambda x: (x['ref_name'], x['pos']))
        
        print(f"✅ Extracted {len(self.tetM_pileup_data)} positions of pileup data")
        
        # Show summary by reference
        ref_counts = {}
        for entry in self.tetM_pileup_data:
            ref_name = entry['ref_name']
            if ref_name not in ref_counts:
                ref_counts[ref_name] = {'count': 0, 'min_pos': float('inf'), 'max_pos': 0}
            ref_counts[ref_name]['count'] += 1
            ref_counts[ref_name]['min_pos'] = min(ref_counts[ref_name]['min_pos'], entry['pos'])
            ref_counts[ref_name]['max_pos'] = max(ref_counts[ref_name]['max_pos'], entry['pos'])
        
        for ref_name, info in ref_counts.items():
            print(f"  {ref_name}: {info['count']} positions ({info['min_pos']} - {info['max_pos']})")
        
        return True
    
    def analyze_with_srst2_algorithm(self):
        """Apply exact SRST2 truncation detection algorithm to tetM data"""
        if not self.tetM_pileup_data:
            print("❌ No tetM pileup data available")
            return False
        
        print(f"\n{'='*80}")
        print("SRST2 TRUNCATION DETECTION ANALYSIS FOR tetM")
        print(f"{'='*80}")
        
        # Group data by reference (in case we have multiple tetM references)
        refs_data = {}
        for entry in self.tetM_pileup_data:
            ref_name = entry['ref_name']
            if ref_name not in refs_data:
                refs_data[ref_name] = []
            refs_data[ref_name].append(entry)
        
        print(f"Analyzing {len(refs_data)} tetM reference(s)...")
        
        all_truncations = {}
        
        for ref_name, positions_data in refs_data.items():
            print(f"\n{'='*60}")
            print(f"ANALYZING REFERENCE: {ref_name}")
            print(f"{'='*60}")
            
            L = self.tetM_info['length']  # Reference length from SRST2 results
            
            print(f"Parameters:")
            print(f"  HIGH_END_FRAC_THRESHOLD: {self.HIGH_END_FRAC_THRESHOLD} (75%)")
            print(f"  VERY_HIGH_END_FRAC: {self.VERY_HIGH_END_FRAC} (90%)")
            print(f"  COVERAGE_DROP_RATIO: {self.COVERAGE_DROP_RATIO} (3x)")
            print(f"  MIN_DEPTH: {self.MIN_DEPTH}")
            print(f"  EDGE_PAD: {self.EDGE_PAD}")
            print(f"  WINDOW_SIZE: {self.WINDOW_SIZE}")
            print(f"  Reference length: {L} bp")
            
            # Track all positions that meet basic criteria
            candidate_positions = []
            
            print(f"\nSTEP-BY-STEP ANALYSIS:")
            
            for i, pos_data in enumerate(positions_data):
                pos = pos_data['pos']
                depth = pos_data['depth']
                bases = pos_data['bases']
                
                # Step 1: Check basic criteria
                meets_min_depth = depth >= self.MIN_DEPTH
                within_edges = not (pos <= self.EDGE_PAD or pos >= (L - self.EDGE_PAD + 1))
                
                if not meets_min_depth or not within_edges:
                    continue
                
                # Count reads ending at this position
                end_count = bases.count('$')
                start_count = bases.count('^')
                frac_end = end_count / float(depth) if depth > 0 else 0.0
                
                # Only show positions with some read ends
                if end_count > 0:
                    print(f"\nPosition {pos}:")
                    print(f"  Depth: {depth}")
                    print(f"  Reads ending ($): {end_count}")
                    print(f"  Reads starting (^): {start_count}")
                    print(f"  End fraction: {frac_end:.3f} ({frac_end*100:.1f}%)")
                    
                    # Check if this position meets the high end fraction threshold
                    meets_threshold = frac_end >= self.HIGH_END_FRAC_THRESHOLD
                    print(f"  ≥75% threshold: {'✅ YES' if meets_threshold else '❌ NO'} ({frac_end:.3f} vs {self.HIGH_END_FRAC_THRESHOLD})")
                    
                    if meets_threshold:
                        candidate_positions.append({
                            'pos': pos,
                            'depth': depth,
                            'end_count': end_count,
                            'frac_end': frac_end,
                            'index': i
                        })
            
            # Analyze candidates for this reference
            truncations_detected = []
            
            if not candidate_positions:
                print(f"\n❌ NO POSITIONS MEET 75% END FRACTION THRESHOLD for {ref_name}")
            else:
                print(f"\n✅ Found {len(candidate_positions)} candidate position(s) for {ref_name}")
                
                for candidate in candidate_positions:
                    pos = candidate['pos']
                    depth = candidate['depth']
                    end_count = candidate['end_count']
                    frac_end = candidate['frac_end']
                    i = candidate['index']
                    
                    print(f"\n--- Analyzing Position {pos} ---")
                    print(f"  End fraction: {frac_end:.3f} ({frac_end*100:.1f}%)")
                    
                    is_likely_truncation = False
                    analysis_reason = ""
                    
                    if frac_end >= self.VERY_HIGH_END_FRAC:
                        # Very high fraction (≥90%) - accept without coverage drop check
                        is_likely_truncation = True
                        analysis_reason = f"Very high end fraction (≥90%): {frac_end:.3f}"
                        print(f"  ✅ TRUNCATION: {analysis_reason}")
                        
                    else:
                        # Moderate fraction (75-90%) - require coverage drop confirmation
                        print(f"  Moderate end fraction (75-90%): {frac_end:.3f}")
                        print(f"  → Checking for coverage drop confirmation...")
                        
                        coverage_drop_confirmed = False
                        drop_details = []
                        
                        # Look at next positions for coverage drop
                        for j in range(i + 1, min(i + 1 + self.WINDOW_SIZE, len(positions_data))):
                            next_pos_data = positions_data[j]
                            next_pos = next_pos_data['pos']
                            next_depth = next_pos_data['depth']
                            
                            # Only consider positions within window
                            if next_pos > pos + self.WINDOW_SIZE:
                                print(f"    Position {next_pos} beyond window (>{pos + self.WINDOW_SIZE})")
                                break
                            
                            # Skip very low coverage positions
                            if next_depth < 5:
                                print(f"    Position {next_pos}: depth {next_depth} too low (<5), skipping")
                                continue
                            
                            # Check for significant coverage drop
                            drop_ratio = depth / next_depth if next_depth > 0 else float('inf')
                            meets_drop_threshold = depth >= next_depth * self.COVERAGE_DROP_RATIO
                            
                            print(f"    Position {next_pos}: depth {next_depth}, ratio {drop_ratio:.1f}x")
                            
                            if meets_drop_threshold:
                                coverage_drop_confirmed = True
                                analysis_reason = f"Coverage drop confirmed: {depth}→{next_depth} ({drop_ratio:.1f}x ≥ {self.COVERAGE_DROP_RATIO}x)"
                                print(f"    ✅ Coverage drop confirmed: {drop_ratio:.1f}x ≥ {self.COVERAGE_DROP_RATIO}x")
                                break
                            else:
                                drop_details.append(f"{next_pos}:{drop_ratio:.1f}x")
                        
                        is_likely_truncation = coverage_drop_confirmed
                        
                        if not coverage_drop_confirmed:
                            analysis_reason = f"Coverage drop NOT confirmed (checked: {', '.join(drop_details)})"
                            print(f"  ❌ NO TRUNCATION: {analysis_reason}")
                    
                    if is_likely_truncation:
                        truncations_detected.append({
                            'ref_name': ref_name,
                            'pos': pos,
                            'frac': frac_end,
                            'depth': depth,
                            'end_count': end_count,
                            'reason': analysis_reason
                        })
            
            # Store results for this reference
            if truncations_detected:
                best_truncation = max(truncations_detected, 
                                    key=lambda x: (x['frac'], x['depth']))
                all_truncations[ref_name] = best_truncation
        
        # Final summary
        print(f"\n{'='*60}")
        print("FINAL RESULTS")
        print(f"{'='*60}")
        
        if all_truncations:
            print(f"✅ TRUNCATION(S) DETECTED in {len(all_truncations)} reference(s)")
            for ref_name, trunc_info in all_truncations.items():
                print(f"\n  Reference: {ref_name}")
                print(f"  Position: {trunc_info['pos']}")
                print(f"  End fraction: {trunc_info['frac']:.3f} ({trunc_info['frac']*100:.1f}%)")
                print(f"  Depth: {trunc_info['depth']}")
                print(f"  Reads ending: {trunc_info['end_count']}")
                print(f"  Reason: {trunc_info['reason']}")
            
            # Return the best overall truncation
            best_overall = max(all_truncations.values(), 
                             key=lambda x: (x['frac'], x['depth']))
            return {
                'truncated': True,
                'reference': best_overall['ref_name'],
                'position': best_overall['pos'],
                'fraction': best_overall['frac'],
                'depth': best_overall['depth'],
                'end_count': best_overall['end_count'],
                'reason': best_overall['reason']
            }
        else:
            print(f"❌ NO TRUNCATION DETECTED")
            print(f"  Reason: No positions passed both 75% threshold AND coverage drop confirmation")
            
            return {
                'truncated': False,
                'reason': 'No positions with ≥75% read ends showed required coverage drop'
            }
    
    def detailed_position_report(self):
        """Generate detailed report of all positions"""
        if not self.tetM_pileup_data:
            return
        
        print(f"\n{'='*80}")
        print("DETAILED POSITION-BY-POSITION REPORT")
        print(f"{'='*80}")
        
        # Group by reference
        refs_data = {}
        for entry in self.tetM_pileup_data:
            ref_name = entry['ref_name']
            if ref_name not in refs_data:
                refs_data[ref_name] = []
            refs_data[ref_name].append(entry)
        
        for ref_name, positions_data in refs_data.items():
            print(f"\nReference: {ref_name}")
            print(f"{'Pos':<6} {'Depth':<6} {'Ends':<6} {'Starts':<7} {'End%':<8} {'≥Min?':<7} {'≥75%?':<7} {'≥90%?':<7} {'Edge?':<7}")
            print("-" * 70)
            
            L = self.tetM_info['length']
            
            for pos_data in positions_data:
                pos = pos_data['pos']
                depth = pos_data['depth']
                bases = pos_data['bases']
                
                end_count = bases.count('$')
                start_count = bases.count('^')
                frac_end = end_count / float(depth) if depth > 0 else 0.0
                
                meets_min_depth = depth >= self.MIN_DEPTH
                meets_75_threshold = frac_end >= self.HIGH_END_FRAC_THRESHOLD
                meets_90_threshold = frac_end >= self.VERY_HIGH_END_FRAC
                within_edges = not (pos <= self.EDGE_PAD or pos >= (L - self.EDGE_PAD + 1))
                
                print(f"{pos:<6} {depth:<6} {end_count:<6} {start_count:<7} {frac_end*100:<7.1f} "
                      f"{'✅' if meets_min_depth else '❌':<7} "
                      f"{'✅' if meets_75_threshold else '❌':<7} "
                      f"{'✅' if meets_90_threshold else '❌':<7} "
                      f"{'✅' if within_edges else '❌':<7}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze tetM truncation using SRST2 algorithm'
    )
    parser.add_argument('results_file', help='SRST2 fullgenes results file')
    parser.add_argument('pileup_file', help='SRST2 pileup file')
    parser.add_argument('--detailed', action='store_true',
                       help='Show detailed position-by-position report')
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = TetMTruncationAnalyzer(args.results_file, args.pileup_file)
    
    # Step 1: Extract tetM from results
    if not analyzer.extract_tetM_from_results():
        return 1
    
    # Step 2: Extract tetM pileup data
    if not analyzer.extract_tetM_pileup_data():
        return 1
    
    # Step 3: Run SRST2 truncation analysis
    result = analyzer.analyze_with_srst2_algorithm()
    
    # Step 4: Detailed report if requested
    if args.detailed:
        analyzer.detailed_position_report()
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    
    if result and result['truncated']:
        print(f"🔴 tetM gene appears to be TRUNCATED")
        print(f"   Reference: {result['reference']}")
        print(f"   Position: {result['position']}")
        print(f"   {result['end_count']} out of {result['depth']} reads end at this position ({result['fraction']*100:.1f}%)")
        print(f"   Algorithm decision: {result['reason']}")
    else:
        print(f"🟢 tetM gene appears to be COMPLETE")
        if result:
            print(f"   Algorithm decision: {result['reason']}")

if __name__ == '__main__':
    main()