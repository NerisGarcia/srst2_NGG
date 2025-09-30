#!/usr/bin/env python3
"""
Script to transform FASTA headers from custom format to AMRFinderPlus format

Input format: >tet(M)_45|GCF_050484505.1.NZ_CP174873.1:323694-325610|tet(M)|NCBIcg_amrfplus
Output format: >0|GCF_050484505.1|NZ_CP174873|1|1|tet(M)|tet(M)-45|tet(M)tetracycline_resistance_ribosomal_protection_protein_Tet(M) NZ_CP174873.1:323694-325610

AMRFinderPlus format fields (separated by |):
1. Protein GI (0 if no GI)
2. Protein accession  
3. DNA accession
4. Fusion gene part number (1 or 2 for fusion, 1 if not fusion)
5. Total fusion parts (2 for fusion, 1 if not fusion)
6. node_id or allele symbol
7. Parent node_id (same as field 5 if not an allele)
8. Resistance mechanism type
9. Protein name with spaces replaced by underscores
"""

import sys
import re
import argparse

def load_amrfinder_fam_tab(fam_tab_file):
    """
    Load AMRFinderPlus fam.tab file and return dictionaries mapping node_id to gene_symbol and family_name.
    """
    gene_symbol_dict = {}
    fam_dict = {}
    try:
        with open(fam_tab_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 18:  # Make sure we have enough columns
                    node_id = parts[0]
                    gene_symbol = parts[2]  # gene_symbol is the 3rd column (index 2)
                    family_name = parts[17]  # family_name is the 18th column (index 17)
                    gene_symbol_dict[node_id] = gene_symbol
                    fam_dict[node_id] = family_name
                    
    except IOError as e:
        print(f"Error reading fam.tab file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return gene_symbol_dict, fam_dict

def transform_header(header_line, node_id_override=None, resistance_mech_override=None, gene_symbol_dict=None, fam_dict=None):
    """Transform a single FASTA header from input format to AMRFinderPlus format"""
    
    # Remove the '>' character
    header = header_line.strip().lstrip('>')
    
    # Split by pipe character
    parts = header.split('|')
    
    if len(parts) < 4:
        print(f"Warning: Header has fewer than 4 parts: {header_line}", file=sys.stderr)
        return header_line
    
    # Extract components from input format
    gene_allele = parts[0]  # e.g., "tet(M)_45"
    location_info = parts[1]  # e.g., "GCF_050484505.1.NZ_CP174873.1:323694-325610"
    gene_name = parts[2]  # e.g., "tet(M)"
    source = parts[3]  # e.g., "NCBIcg_amrfplus"
    
    # Parse gene and allele number
    if '_' in gene_allele:
        gene_base, allele_num = gene_allele.rsplit('_', 1)
    else:
        gene_base = gene_allele
        allele_num = "0"
    
    # Parse location information - handle both simple and complex formats
    # Simple format: NG_048253.1:101-2020
    # Complex format: GCF_900639655.1.NZ_LR135435.1:2426014-2427930
    if ':' in location_info:
        accession_part, coords = location_info.split(':', 1)
        
        # Check if it's the complex format with GCF
        if '.NZ_' in accession_part:
            gcf_acc, dna_acc = accession_part.split('.NZ_', 1)
            dna_acc = 'NZ_' + dna_acc
        elif '.NC_' in accession_part:
            gcf_acc, dna_acc = accession_part.split('.NC_', 1)
            dna_acc = 'NC_' + dna_acc
        else:
            # Simple format - use the same for both
            gcf_acc = accession_part
            dna_acc = accession_part
    else:
        # No coordinates
        gcf_acc = location_info
        dna_acc = location_info
        coords = ""
    
    # Build AMRFinderPlus format
    # Fields: protein_gi|protein_acc|dna_acc|fusion_part|fusion_total|node_id|parent_node|mechanism|protein_name coordinates
    
    protein_gi = "0"  # No GI provided
    protein_acc = gcf_acc  # Use GCF as protein accession
    dna_accession = dna_acc
    fusion_part = "1"
    fusion_total = "1" 
    
    # Determine node_id and gene_symbol
    if node_id_override:
        node_id = node_id_override
    else:
        node_id = gene_name  # e.g., "tet(M)"
    
    # Get gene_symbol from fam.tab
    if gene_symbol_dict and node_id in gene_symbol_dict:
        gene_symbol = gene_symbol_dict[node_id]
    else:
        gene_symbol = node_id  # fallback to node_id if not found
    
    parent_node = f"{gene_name}-{allele_num}"  # e.g., "tet(M)-45"
    
    # Determine resistance mechanism
    if resistance_mech_override:
        mechanism_type = resistance_mech_override
    elif fam_dict and node_id in fam_dict:
        mechanism_type = fam_dict[node_id].replace(' ', '_')
    else:
        # Fallback - use gene name
        mechanism_type = f"{gene_name}_antimicrobial_resistance_protein"
        print(f"Warning: No mechanism found for {node_id}, using fallback", file=sys.stderr)
    
    # Reconstruct coordinates
    coordinates = f"{dna_acc}:{coords}" if coords else dna_acc
    
    # Build the final header - Changed order: parent_node before gene_symbol
    amr_header = f">{protein_gi}|{protein_acc}|{dna_accession}|{fusion_part}|{fusion_total}|{parent_node}|{gene_symbol}|{mechanism_type}"
    
    if coordinates:
        amr_header += f" {coordinates}"
    
    return amr_header

def process_fasta_file(input_file, output_file, node_id_override=None, resistance_mech_override=None, gene_symbol_dict=None, fam_dict=None):
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
                        new_header = transform_header(line, node_id_override, resistance_mech_override, gene_symbol_dict, fam_dict)
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
    parser = argparse.ArgumentParser(description='Transform FASTA headers to AMRFinderPlus format')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    
    # Manual override options
    parser.add_argument('--node-id', help='Manual node_id override')
    parser.add_argument('--resistance-mech', help='Manual resistance mechanism override')
    
    # AMRFinderPlus fam.tab file option
    parser.add_argument('--amrfplus-fam-tab', help='Path to AMRFinderPlus fam.tab file for automatic lookup')
    
    parser.add_argument('--test', action='store_true', help='Test with example header')
    
    args = parser.parse_args()
    
    # Load fam.tab file if provided
    gene_symbol_dict = None
    fam_dict = None
    if args.amrfplus_fam_tab:
        gene_symbol_dict, fam_dict = load_amrfinder_fam_tab(args.amrfplus_fam_tab)
        print(f"Loaded {len(fam_dict)} entries from fam.tab file", file=sys.stderr)
    
    if args.test:
        # Test with the complex format (original example)
        test_header1 = ">tet(M)_45|GCF_050484505.1.NZ_CP174873.1:323694-325610|tet(M)|NCBIcg_amrfplus"
        result1 = transform_header(test_header1, args.node_id, args.resistance_mech, gene_symbol_dict, fam_dict)
        print("Input (complex):", test_header1)
        print("Output:", result1)
        
        # Test with the simple format (after sed preprocessing)
        test_header2 = ">tet(M)_1|NG_048253.1:101-2020|tet(M)|AMR_ref_gene_catalog"
        result2 = transform_header(test_header2, args.node_id, args.resistance_mech, gene_symbol_dict, fam_dict)
        print("\nInput (simple):", test_header2)
        print("Output:", result2)
        
        # Test with another complex format
        test_header3 = ">tet(M)_43|GCF_900639655.1.NZ_LR135435.1:2426014-2427930|tet(M)|NCBIcg_amrfplus"
        result3 = transform_header(test_header3, args.node_id, args.resistance_mech, gene_symbol_dict, fam_dict)
        print("\nInput (complex 2):", test_header3)
        print("Output:", result3)
        return
    
    process_fasta_file(args.input_file, args.output, args.node_id, args.resistance_mech, gene_symbol_dict, fam_dict)

if __name__ == "__main__":
    main()
