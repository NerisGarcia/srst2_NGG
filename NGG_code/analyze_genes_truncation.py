#!/usr/bin/env python3

import argparse
import os
import sys
from collections import defaultdict
import pandas as pd

class GenesAnalyzer:
    def __init__(self, results_file, pileup_file):
        self.results_file = results_file
        self.pileup_file = pileup_file
        
        # SRST2 parameters
        self.HIGH_END_FRAC_THRESHOLD = 0.75   
        self.VERY_HIGH_END_FRAC = 0.90        
        self.COVERAGE_DROP_RATIO = 3.0        
        self.EDGE_PAD = 10                    
        self.WINDOW_SIZE = 10                 
        self.LOW_COV_THRESHOLD = 0.4  # 40% of mean depth
        
        # Store data
        self.pileup_data = defaultdict(list)
        self.ref_to_gene = {}  # Map pileup references to gene names
        self.mean_depth = None  # Will store mean depth of non-truncated genes
        
    def analyze_genes(self):
        """Analyze all genes and their variants for truncations"""
        print("Reading SRST2 results...")
        
        # Read results file to get all genes and variants
        genes_info = []
        with open(self.results_file) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= len(header):
                    gene_info = dict(zip(header, parts))
                    genes_info.append(gene_info)
        
        # Calculate mean depth from all genes
        depths = []
        for gene_info in genes_info:
            if 'depth' in gene_info:
                try:
                    depth = float(gene_info['depth'])
                    depths.append(depth)
                except ValueError:
                    continue
        
        if depths:
            self.mean_depth = sum(depths) / len(depths)
            print(f"Mean depth across genes: {self.mean_depth:.2f}x")
            print(f"Low coverage threshold (<{self.LOW_COV_THRESHOLD*100}%): {self.mean_depth * self.LOW_COV_THRESHOLD:.2f}x")
        
        # Read pileup data
        print("Reading pileup data...")
        self._read_pileup_data()
        
        # Analyze each gene
        results = []
        for gene_info in genes_info:
            gene = gene_info['gene']
            allele = gene_info['allele']
            gene_depth = float(gene_info['depth'])
            
            # Find matching references in pileup
            matching_refs = self._find_matching_refs(gene, allele)
            
            if matching_refs:
                # Analyze each reference
                for ref in matching_refs:
                    truncation_data = self._check_truncation(ref, int(gene_info['length']))
                    
                    results.append({
                        'Gene': gene,
                        'Allele': allele,
                        'Reference': ref,
                        'Length': int(gene_info['length']),
                        'Coverage': float(gene_info['coverage']),
                        'Depth': gene_depth,
                        'Low_Coverage': gene_depth < (self.mean_depth * self.LOW_COV_THRESHOLD) if self.mean_depth else False,
                        'Truncated': truncation_data['truncated'],
                        'Position': truncation_data.get('position', '-'),
                        'End_Fraction': f"{truncation_data.get('fraction', 0)*100:.1f}%" if truncation_data.get('fraction') else '-',
                        'Depth_at_pos': truncation_data.get('depth', '-'),
                        'Clipped_reads': truncation_data.get('end_count', '-'),
                        'Reason': truncation_data.get('reason', '-')
                    })
            else:
                print(f"\nWarning: No pileup data found for {gene} ({allele})")
        
        return pd.DataFrame(results)
    
    def _read_pileup_data(self):
        """Read full pileup file and organize by reference"""
        self.pileup_data.clear()
        
        with open(self.pileup_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    ref = parts[0]
                    self.pileup_data[ref].append({
                        'pos': int(parts[1]),
                        'ref_base': parts[2],
                        'depth': int(parts[3]),
                        'bases': parts[4]
                    })
    
    def _find_matching_refs(self, gene, allele):
        """Find pileup references matching gene/allele"""
        matching_refs = []
        
        # Convert gene and allele names to match pileup format
        gene_patterns = [
            gene.lower(),
            gene.replace('(', '').replace(')', '').lower(),
            allele.lower(),
            allele.replace('(', '').replace(')', '').lower()
        ]
        
        for ref in self.pileup_data.keys():
            ref_lower = ref.lower()
            if any(pattern in ref_lower for pattern in gene_patterns):
                matching_refs.append(ref)
                self.ref_to_gene[ref] = gene
        
        return matching_refs
    
    def _check_truncation(self, ref, length):
        """Check if a gene reference shows truncation"""
        if ref not in self.pileup_data:
            return {'truncated': False, 'reason': 'No pileup data found'}
        
        positions = self.pileup_data[ref]
        if not positions:
            return {'truncated': False, 'reason': 'No positions in pileup'}
        
        # Look for truncation positions
        truncation_candidates = []
        
        for i, pos_data in enumerate(positions):
            pos = pos_data['pos']
            depth = pos_data['depth']
            bases = pos_data['bases']
            
            # Skip only edge positions
            if pos <= self.EDGE_PAD or pos >= (length - self.EDGE_PAD + 1):
                continue
            
            # Count reads ending at this position ($)
            end_count = bases.count('$')
            frac_end = end_count / float(depth) if depth > 0 else 0.0
            
            if frac_end >= self.HIGH_END_FRAC_THRESHOLD:
                if frac_end >= self.VERY_HIGH_END_FRAC:
                    truncation_candidates.append({
                        'pos': pos,
                        'depth': depth,
                        'end_count': end_count,
                        'fraction': frac_end,
                        'reason': f"≥90% reads end ({frac_end*100:.1f}%)"
                    })
                else:
                    # Check coverage drop
                    max_window = min(i + self.WINDOW_SIZE, len(positions))
                    for j in range(i + 1, max_window):
                        next_pos = positions[j]['pos']
                        next_depth = positions[j]['depth']
                        
                        if next_pos > pos + self.WINDOW_SIZE:
                            break
                            
                        if next_depth > 0:
                            drop_ratio = depth / next_depth
                            
                            if drop_ratio >= self.COVERAGE_DROP_RATIO:
                                truncation_candidates.append({
                                    'pos': pos,
                                    'depth': depth,
                                    'end_count': end_count,
                                    'fraction': frac_end,
                                    'reason': f"75-90% ends + {drop_ratio:.1f}x coverage drop"
                                })
                                break
        
        if truncation_candidates:
            # Select best candidate (highest fraction, then highest depth)
            best = max(truncation_candidates, key=lambda x: (x['fraction'], x['depth']))
            return {
                'truncated': True,
                'position': best['pos'],
                'depth': best['depth'],
                'end_count': best['end_count'],
                'fraction': best['fraction'],
                'reason': best['reason']
            }
        
        return {'truncated': False, 'reason': 'No truncation detected'}

def main():
    parser = argparse.ArgumentParser(description='Analyze genes for truncations')
    parser.add_argument('results_file', help='SRST2 fullgenes results file')
    parser.add_argument('pileup_file', help='SRST2 pileup file')
    parser.add_argument('--output', help='Output CSV file')
    
    args = parser.parse_args()
    
    analyzer = GenesAnalyzer(args.results_file, args.pileup_file)
    results_df = analyzer.analyze_genes()
    
    # Print results
    print("\nResults:")
    print(results_df.to_string(index=False))
    
    # Save to CSV if requested
    if args.output:
        results_df.to_csv(args.output, index=False)
        print(f"\nResults saved to: {args.output}")

if __name__ == '__main__':
    main()