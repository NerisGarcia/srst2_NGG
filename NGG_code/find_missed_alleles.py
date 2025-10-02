#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from collections import defaultdict
import logging

def parse_srst2_fullgenes(fullgenes_file):
    """Parse SRST2 fullgenes report to get reported alleles and their uncertainty"""
    reported = {}
    
    with open(fullgenes_file) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            gene = fields[2]
            allele = fields[3]
            uncertainty = fields[7]  # uncertainty column
            diffs = fields[6]      # differences column
            
            reported[allele] = {
                'gene': gene,
                'allele': allele,
                'has_uncertainty': bool(uncertainty.strip()),
                'has_diffs': bool(diffs.strip())
            }
    
    return reported

def analyze_gene_coverage(pileup_file, min_coverage=60):
    """
    Analyze pileup for allele detection using SRST2-like criteria:
    1. Coverage (% positions with depth > 0)
    2. Depth compared to mean
    3. Truncation detection
    4. Reports divergence for information
    5. Low coverage detection
    """
    # Get reference lengths
    ref_lengths = {}
    with open(pileup_file) as f:
        for line in f:
            ref = line.strip().split('\t')[0]
            if ref not in ref_lengths:
                ref_lengths[ref] = 0
            ref_lengths[ref] += 1
    
    # Initialize stats tracking
    allele_stats = defaultdict(lambda: {
        'depths': [],
        'positions_with_depth': 0,
        'total_pos': 0,
        'mismatches': 0,
        'end_counts': defaultdict(int)  # For truncation detection
    })
    
    # First pass: collect depth statistics
    total_depth = 0
    total_positions = 0
    
    with open(pileup_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            allele = fields[0]  # The allele name is directly in the first column
            pos = int(fields[1])
            depth = int(fields[3])
            
            total_depth += depth
            total_positions += 1
            
            allele_stats[allele]['depths'].append(depth)
            allele_stats[allele]['total_pos'] += 1
            if depth > 0:
                allele_stats[allele]['positions_with_depth'] += 1
            
            # Count mismatches and read ends
            if len(fields) > 4:
                bases = fields[4].upper()
                ref_base = fields[2].upper()
                allele_stats[allele]['mismatches'] += sum(1 for b in bases if b not in ['.', ',', ref_base])
                allele_stats[allele]['end_counts'][pos] = bases.count('$')
    
    # Calculate overall mean depth for low coverage detection
    mean_depth = float(total_depth)/total_positions if total_positions > 0 else 0
    low_cov_threshold = mean_depth * 0.1  # 10% of mean depth
    
    # Process final results
    results = {}
    for allele_id, stats in allele_stats.items():
        if stats['total_pos'] > 0:
            # Calculate basic metrics
            coverage = 100.0 * stats['positions_with_depth'] / stats['total_pos']
            mean_allele_depth = sum(stats['depths']) / float(len(stats['depths']))
            divergence = 100.0 * stats['mismatches'] / stats['positions_with_depth'] if stats['positions_with_depth'] > 0 else 0
            
            # Detect truncation
            is_truncated = False
            max_end_frac = 0
            truncation_pos = None
            
            for pos, end_count in stats['end_counts'].items():
                if pos <= 10 or pos >= stats['total_pos'] - 10:  # Skip edges
                    continue
                depth = stats['depths'][pos-1]  # pos is 1-based
                if depth > 0:
                    end_frac = float(end_count)/depth
                    if end_frac >= 0.9:  # Very high threshold
                        is_truncated = True
                        max_end_frac = end_frac
                        truncation_pos = pos
                        break
                    elif end_frac >= 0.7:  # Check for coverage drop
                        next_depths = stats['depths'][pos:pos+10]
                        if next_depths:
                            avg_next = sum(next_depths)/len(next_depths)
                            if depth/avg_next >= 3.0:
                                is_truncated = True
                                max_end_frac = end_frac
                                truncation_pos = pos
                                break
            
            # Store results if coverage passes minimum
            if coverage >= min_coverage:
                results[allele] = {
                    'coverage': coverage,
                    'mean_depth': mean_allele_depth,
                    'divergence': divergence,
                    'is_truncated': is_truncated,
                    'truncation_pos': truncation_pos,
                    'truncation_frac': max_end_frac if is_truncated else 0,
                    'is_low_depth': mean_allele_depth < low_cov_threshold,
                    'depth_vs_mean': mean_allele_depth/mean_depth if mean_depth > 0 else 0
                }
    
    return results

def get_gene_from_allele(allele_name):
    """Extract gene name from allele name (between underscores)"""
    parts = allele_name.split('_')
    if len(parts) >= 3:  # Must have at least one underscore-separated part
        return '_'.join(parts[:-1])  # Return everything except the last part
    return allele_name

def main():
    parser = argparse.ArgumentParser(description='Find potentially missed alleles in SRST2 results')
    parser.add_argument('--fullgenes', required=True, help='SRST2 fullgenes report')
    parser.add_argument('--pileup', required=True, help='SRST2 pileup file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--min_coverage', type=float, default=60,
                       help='Minimum coverage to consider allele present (default: 60)')
    args = parser.parse_args()
    
    # Get reported alleles and their flags
    reported_alleles = parse_srst2_fullgenes(args.fullgenes)

    # Analyze pileup
    allele_results = analyze_gene_coverage(args.pileup, args.min_coverage)
    
    # Write results
    with open(args.output, 'w') as out:
        out.write("\t".join([
            "Sample", "Gene", "Allele", "Status", "Coverage", "Mean_Depth", 
            "Depth_vs_Mean", "Divergence", "Is_Truncated", "Truncation_Pos",
            "Is_Low_Depth", "SRST2_Status", "Has_Uncertainty", "Has_Differences"
        ]) + "\n")
        
        sample = os.path.basename(args.pileup).split('.')[0]
        
        for allele, stats in allele_results.items():
            gene = get_gene_from_allele(allele)
            
            # Determine status and reporting flags
            if allele in reported_alleles:
                reported = reported_alleles[allele]
                if not reported['has_uncertainty'] and not reported['has_diffs']:
                    status = "allele_already_detected"
                    srst2_status = "Perfect_Match"
                else:
                    status = "allele_reported_with_flags"
                    srst2_status = "With_Flags"
                has_uncertainty = str(reported['has_uncertainty'])
                has_diffs = str(reported['has_diffs'])
            else:
                status = "new_allele_detected"
                srst2_status = "Not_Reported"
                has_uncertainty = "-"
                has_diffs = "-"
            
            # Write result with all details
            out.write("\t".join([
                sample,
                gene,
                allele,
                status,
                f"{stats['coverage']:.1f}",
                f"{stats['mean_depth']:.1f}",
                f"{stats['depth_vs_mean']:.2f}",
                f"{stats['divergence']:.1f}",
                str(stats['is_truncated']),
                str(stats['truncation_pos'] if stats['is_truncated'] else "-"),
                str(stats['is_low_depth']),
                srst2_status,
                has_uncertainty,
                has_diffs
            ]) + "\n")

if __name__ == '__main__':
    main()