#!/usr/bin/env python3
"""
Script to simplify AMRFinderPlus FASTA headers for SRST2 clustering

Input format: >99030340|ABF61437.1|DQ517331.2|1|1|ncrY|ncrY|nickel_resistance_OB_fold_protein_NcrY DQ517331.2:3054-3476
Output format: >tet(M)_tet(M)-71 ABF61437.1_DQ517331.2:3054-3476

This creates simplified headers suitable for CD-HIT clustering, keeping only:
- node_id (gene name)
- allele number  
- protein accession
- location (DNA accession + coordinates)

The format is designed to work with SRST2's gene-allele parsing which splits on "-"
"""

import sys
import re
import argparse

def parse_allele_info(parent_node, node_id):
    """Extract allele number from parent_node"""
    if '-' in parent_node and parent_node != node_id:
        # Format like "tet(M)-1"
        return parent_node.split('-')[-1]
    else:
        return "1"

def transform_header_for_clustering(header_line):
    """Transform AMRFinderPlus header to simplified format for clustering"""
    
    # Remove the '>' character
    header = header_line.strip().lstrip('>')
    
    # Split by pipe character first
    if '|' not in header:
        print(f"Warning: Header doesn't contain pipe separators: {header_line}", file=sys.stderr)
        return header_line
    
    # Split header into pipe-separated parts and additional info (coordinates)
    if ' ' in header:
        pipe_part, coordinates = header.split(' ', 1)
    else:
        pipe_part = header
        coordinates = ""
    
    parts = pipe_part.split('|')
    
    if len(parts) < 7:
        print(f"Warning: Header has fewer than 7 pipe-separated parts: {header_line}", file=sys.stderr)
        return header_line
    
    # Extract components - Updated for new field order
    protein_gi = parts[0]
    protein_acc = parts[1] 
    dna_acc = parts[2]
    fusion_part = parts[3]
    fusion_total = parts[4]
    parent_node = parts[5]  # allele symbol (e.g., "tet(M)-1")
    gene_symbol = parts[6]  # gene symbol (e.g., "tet(M)")
    
    # Get allele number from parent_node
    allele_number = parse_allele_info(parent_node, gene_symbol)
    
    # Build simplified header: >gene_symbol_allele additional_info
    # This format ensures SRST2 can properly parse gene-allele using the first "-" in the allele part
    
    # Gene name (gene_symbol)
    gene_name = gene_symbol
    
    # Full allele name (parent_node or constructed)
    if parent_node != gene_symbol:
        allele_name = parent_node  # e.g., "tet(M)-71"
    else:
        allele_name = f"{gene_symbol}-{allele_number}"  # e.g., "ncrY-1"
    
    # Main identifier: gene_symbol_allele
    main_id = f"{gene_name}_{allele_name}"
    
    # Additional info components
    additional_info = []
    
    # Add protein accession if available
    if protein_acc and protein_acc != "0":
        additional_info.append(protein_acc)
    
    # Add location (DNA accession + coordinates)
    if coordinates:
        location = coordinates
    else:
        location = dna_acc if dna_acc else ""
    
    if location:
        additional_info.append(location)
    
    # Build final header
    simplified_header = f">{main_id}"
    if additional_info:
        simplified_header += " " + "_".join(additional_info)
    
    return simplified_header

def process_fasta_file(input_file, output_file):
    """Process entire FASTA file, transforming headers"""
    
    try:
        with open(input_file, 'r') as infile:
            if output_file:
                outfile = open(output_file, 'w')
            else:
                outfile = sys.stdout
            
            try:
                for line in infile:
                    line = line.rstrip('\n\r')
                    
                    if line.startswith('>'):
                        # Transform header
                        new_header = transform_header_for_clustering(line)
                        outfile.write(new_header + '\n')
                    else:
                        # Copy sequence line as-is
                        outfile.write(line + '\n')
                        
            finally:
                if output_file:
                    outfile.close()
                    
    except IOError as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Simplify AMRFinderPlus headers for SRST2 clustering')
    parser.add_argument('input_file', nargs='?', help='Input FASTA file with AMRFinderPlus headers')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('--test', action='store_true', help='Test with example headers')
    
    args = parser.parse_args()
    
    if args.test:
        # Test with various example headers - Updated for new format
        test_headers = [
            ">99030340|ABF61437.1|DQ517331.2|1|1|ncrY-1|ncrY|nickel_resistance_OB_fold_protein_NcrY DQ517331.2:3054-3476",
            ">991872253|AMG40538.1|CP014064.2|1|1|sel-1|sel|staphylococcal_enterotoxin_type_L CP014064.2:7864-7142",
            ">0|NG_048253.1|NG_048253.1|1|1|tet(M)-1|tet(M)|tetracycline_resistance_ribosomal_protection_protein_Tet(M) NG_048253.1:101-2020",
            ">0|NG_048252.1|NG_048252.1|1|1|tet(M)-2|tet(M)|tetracycline_resistance_ribosomal_protection_protein_Tet(M) NG_048252.1:101-2020"
        ]
        
        for i, test_header in enumerate(test_headers, 1):
            result = transform_header_for_clustering(test_header)
            print(f"Test {i}:")
            print(f"Input:  {test_header}")
            print(f"Output: {result}")
            print()
        return
    
    if not args.input_file:
        parser.error("input_file is required when not using --test")
    
    process_fasta_file(args.input_file, args.output)

if __name__ == "__main__":
    main()
