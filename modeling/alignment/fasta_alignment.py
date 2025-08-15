#!/usr/bin/env python3
"""
Protein Sequence Alignment and Residue Mapping using BioPython
This script aligns two protein sequences from FASTA files and creates residue mapping dictionaries.
"""

import os
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
from Bio import SeqIO
from Bio import Align
import argparse
import json

matrix = substitution_matrices.load("BLOSUM62")

def read_fasta_sequence(fasta_file):
    """
    Read a single protein sequence from a FASTA file
    Handles multi-line sequences properly
    Returns: (sequence_id, sequence_string)
    """
    try:
        with open(fasta_file, 'r') as f:
            records = list(SeqIO.parse(f, 'fasta'))
            
            if not records:
                raise ValueError(f"No sequences found in {fasta_file}")
            
            if len(records) > 1:
                print(f"Warning: {fasta_file} contains {len(records)} sequences. Using the first one: {records[0].id}")
            
            # Get the first (or only) record
            record = records[0]
            sequence = str(record.seq).upper().strip()
            
            # Remove any remaining whitespace or newlines from sequence
            sequence = ''.join(sequence.split())
            
            # Validate that sequence contains only valid amino acid characters
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY*-')
            invalid_chars = set(sequence) - valid_aa
            if invalid_chars:
                print(f"Warning: Found potentially invalid characters in {fasta_file}: {invalid_chars}")
                # Remove invalid characters
                sequence = ''.join(c for c in sequence if c in valid_aa)
            
            if not sequence:
                raise ValueError(f"Empty sequence in {fasta_file}")
            
            print(f"Successfully read sequence {record.id}: {len(sequence)} residues")
            return record.id, sequence
            
    except FileNotFoundError:
        print(f"Error: File {fasta_file} not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)

def create_residue_mappings(aligned_seq_a, aligned_seq_b):
    """
    Create residue mapping dictionaries between two aligned sequences.
    
    Returns:
        dict: A_to_B mapping {pos_in_A: (residue_A, pos_in_B, residue_B)}
        dict: B_to_A mapping {pos_in_B: (residue_B, pos_in_A, residue_A)}
    """
    A_to_B = {}  # {position_in_A: (residue_A, position_in_B, residue_B)}
    B_to_A = {}  # {position_in_B: (residue_B, position_in_A, residue_A)}
    
    # Track positions in original sequences (1-indexed)
    pos_a = 0  # Current position in original sequence A
    pos_b = 0  # Current position in original sequence B
    
    # Iterate through aligned sequences
    for i in range(len(aligned_seq_a)):
        res_a = aligned_seq_a[i]
        res_b = aligned_seq_b[i]
        
        # Handle sequence A residue
        if res_a != '-':
            pos_a += 1
            
            # Map A to B
            if res_b != '-':
                # Both have residues - direct mapping
                A_to_B[pos_a] = (res_a, pos_b + 1, res_b)
            else:
                # A has residue, B has gap
                A_to_B[pos_a] = (res_a, '-', '-')
        
        # Handle sequence B residue
        if res_b != '-':
            pos_b += 1
            
            # Map B to A
            if res_a != '-':
                # Both have residues - direct mapping
                B_to_A[pos_b] = (res_b, pos_a, res_a)
            else:
                # B has residue, A has gap
                B_to_A[pos_b] = (res_b, '-', '-')
    
    return A_to_B, B_to_A


def perform_alignment_seqs(seq1_id, seq1_seq, seq2_id, seq2_seq, output_prefix='alignment'):
    """
    Perform sequence alignment using BioPython's new Align module
    """
    print(f"Using BioPython Align module for alignment...")
    
    # Create aligner object
    aligner = Align.PairwiseAligner()
    
    # Set alignment parameters
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'  # Global alignment
    
    # Perform alignment
    # alignments = aligner.align(seq1_seq, seq2_seq)
    alignments = pairwise2.align.globalds(seq1_seq, seq2_seq, matrix, -10, -0.5)
    
    # Get the best alignment
    best_alignment = alignments[0]
    aligned_seq_a = best_alignment[0]
    aligned_seq_b = best_alignment[1]
    score = best_alignment.score

    # print(best_alignment[0])
    # print(best_alignment[1])
    
    return aligned_seq_a, aligned_seq_b, score

def perform_alignment(seq1_file, seq2_file, output_prefix='alignment'):
    """
    Perform sequence alignment using BioPython and create mapping dictionaries
    """
    # Read sequences from FASTA files
    seq1_id, seq1_seq = read_fasta_sequence(seq1_file)
    seq2_id, seq2_seq = read_fasta_sequence(seq2_file)
    
    print(f"Protein A: {seq1_id} (length: {len(seq1_seq)})")
    print(f"Protein B: {seq2_id} (length: {len(seq2_seq)})")

    # print(seq1_seq)
    # print(seq2_seq)
    
    try:
        aligned_seq_a, aligned_seq_b, score = perform_alignment_seqs(
            seq1_id, seq1_seq, seq2_id, seq2_seq, output_prefix)
        
        print(f"\nAlignment completed successfully!")
        print(f"Alignment score: {score:.2f}")

        # print(aligned_seq_a)
        # print(aligned_seq_b)
        
        # Write alignment results to files
        output_alignment = f"{output_prefix}_alignment.txt"
        output_fasta = f"{output_prefix}_aligned.fasta"
        
        # Write alignment in readable format
        with open(output_alignment, 'w') as f:
            f.write(f"Alignment between {seq1_id} and {seq2_id}\n")
            f.write(f"Score: {score:.2f}\n\n")
            
            # Split alignment into 100-character chunks
            chunk_size = 100
            for start in range(0, len(aligned_seq_a), chunk_size):
                end = min(start + chunk_size, len(aligned_seq_a))
                
                # Get sequence chunks
                chunk_a = aligned_seq_a[start:end]
                chunk_b = aligned_seq_b[start:end]
                
                # Create match line for this chunk
                match_line = ""
                for i in range(len(chunk_a)):
                    if chunk_a[i] == chunk_b[i]:
                        match_line += "|"
                    elif chunk_a[i] == '-' or chunk_b[i] == '-':
                        match_line += " "
                    else:
                        match_line += "."
                
                # Write this chunk
                f.write(f"A: {chunk_a}\n")
                f.write(f"   {match_line}\n")
                f.write(f"B: {chunk_b}\n")
                f.write(f"\n")  # Blank line between chunks
        
        # Write FASTA format
        with open(output_fasta, 'w') as f:
            f.write(f">{seq1_id}_aligned\n{aligned_seq_a}\n")
            f.write(f">{seq2_id}_aligned\n{aligned_seq_b}\n")
        
        print(f"Output files generated:")
        print(f"  - Alignment text: {output_alignment}")
        print(f"  - Aligned FASTA: {output_fasta}")
        
        # Display partial alignment for verification
        # print(f"\nAligned sequences (first 80 characters):")
        # print(f"A: {aligned_seq_a[:80]}")
        # print(f"B: {aligned_seq_b[:80]}")
        # if len(aligned_seq_a) > 80:
        #     print(f"   ... (total length: {len(aligned_seq_a)})")
        
        # Create residue mapping dictionaries
        A_to_B, B_to_A = create_residue_mappings(aligned_seq_a, aligned_seq_b)
        
        # Save mapping dictionaries
        mapping_file_a_to_b = f"{output_prefix}_A_to_B_mapping.json"
        mapping_file_b_to_a = f"{output_prefix}_B_to_A_mapping.json"
        
        # Convert mappings to JSON-serializable format
        A_to_B_json = {str(k): {'res_A': v[0], 'pos_B': v[1], 'res_B': v[2]} for k, v in A_to_B.items()}
        B_to_A_json = {str(k): {'res_B': v[0], 'pos_A': v[1], 'res_A': v[2]} for k, v in B_to_A.items()}
        
        with open(mapping_file_a_to_b, 'w') as f:
            json.dump(A_to_B_json, f, indent=2)
        
        with open(mapping_file_b_to_a, 'w') as f:
            json.dump(B_to_A_json, f, indent=2)
        
        print(f"  - A to B mapping: {mapping_file_a_to_b}")
        print(f"  - B to A mapping: {mapping_file_b_to_a}")
        
        # Print alignment statistics
        identity = calculate_identity(aligned_seq_a, aligned_seq_b)
        similarity = calculate_similarity(aligned_seq_a, aligned_seq_b)
        gaps_a = aligned_seq_a.count('-')
        gaps_b = aligned_seq_b.count('-')
        
        print(f"\nAlignment statistics:")
        print(f"  - Alignment length: {len(aligned_seq_a)}")
        print(f"  - Sequence identity: {identity:.2f}%")
        print(f"  - Sequence similarity: {similarity:.2f}%")
        print(f"  - Gaps in A: {gaps_a}")
        print(f"  - Gaps in B: {gaps_b}")
        
        return A_to_B, B_to_A, seq1_id, seq2_id
        
    except Exception as e:
        print(f"Error during alignment: {e}")
        sys.exit(1)

def calculate_identity(seq1_aligned, seq2_aligned):
    """
    Calculate sequence identity percentage from aligned sequences
    """
    matches = 0
    total = 0
    
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] != '-' and seq2_aligned[i] != '-':
            total += 1
            if seq1_aligned[i] == seq2_aligned[i]:
                matches += 1
    
    if total > 0:
        return (matches / total) * 100
    else:
        return 0.0

def calculate_similarity(seq1_aligned, seq2_aligned):
    """
    Calculate sequence similarity percentage (including conservative substitutions)
    """
    # Simple similarity - you could enhance this with BLOSUM matrix
    similar_groups = [
        set(['A', 'G']),  # Small
        set(['I', 'L', 'V']),  # Hydrophobic
        set(['F', 'W', 'Y']),  # Aromatic
        set(['D', 'E']),  # Acidic
        set(['K', 'R']),  # Basic
        set(['S', 'T']),  # Polar
        set(['N', 'Q']),  # Amide
    ]
    
    matches = 0
    total = 0
    
    for i in range(len(seq1_aligned)):
        if seq1_aligned[i] != '-' and seq2_aligned[i] != '-':
            total += 1
            if seq1_aligned[i] == seq2_aligned[i]:
                matches += 1
            else:
                # Check if they're in the same similarity group
                for group in similar_groups:
                    if seq1_aligned[i] in group and seq2_aligned[i] in group:
                        matches += 1
                        break
    
    if total > 0:
        return (matches / total) * 100
    else:
        return 0.0

def query_site_mapping(A_to_B, B_to_A, seq1_id, seq2_id, site_query):
    """
    Query specific site mapping
    site_query format: A32 or B15
    """
    print(f"\nSite mapping query: {site_query}")
    
    if site_query.upper().startswith('A'):
        # Query protein A position
        try:
            pos = int(site_query[1:])
            if pos in A_to_B:
                res_a, pos_b, res_b = A_to_B[pos]
                print(f"Protein A position {pos} ({res_a}) maps to:")
                if pos_b != '-':
                    print(f"  Protein B position {pos_b} ({res_b})")
                else:
                    print(f"  Gap in protein B (no corresponding residue)")
            else:
                print(f"Position {pos} not found in protein A (length: {max(A_to_B.keys()) if A_to_B else 0})")
        except ValueError:
            print(f"Invalid position format: {site_query}")
    
    elif site_query.upper().startswith('B'):
        # Query protein B position
        try:
            pos = int(site_query[1:])
            if pos in B_to_A:
                res_b, pos_a, res_a = B_to_A[pos]
                print(f"Protein B position {pos} ({res_b}) maps to:")
                if pos_a != '-':
                    print(f"  Protein A position {pos_a} ({res_a})")
                else:
                    print(f"  Gap in protein A (no corresponding residue)")
            else:
                print(f"Position {pos} not found in protein B (length: {max(B_to_A.keys()) if B_to_A else 0})")
        except ValueError:
            print(f"Invalid position format: {site_query}")
    
    else:
        print(f"Invalid site query format: {site_query}")
        print("Use format: A32 (for protein A position 32) or B15 (for protein B position 15)")

def main():
    parser = argparse.ArgumentParser(description='Align two protein sequences using BioPython and create residue mappings')
    parser.add_argument('-i', '--input', nargs=2, required=True, metavar=('FASTA1', 'FASTA2'),
                       help='Two FASTA files to align')
    parser.add_argument('-s', '--site', 
                       help='Query specific site mapping (e.g., A32 for position 32 in protein A)')
    parser.add_argument('-o', '--output', default='alignment', 
                       help='Output prefix for alignment and mapping files (default: alignment)')
    parser.add_argument('--verbose', action='store_true', 
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    seq1_file, seq2_file = args.input
    
    # Check if input files exist
    for fasta_file in [seq1_file, seq2_file]:
        if not os.path.exists(fasta_file):
            print(f"Error: File {fasta_file} not found")
            sys.exit(1)
    
    print("Starting protein sequence alignment with BioPython...")
    print(f"Input files: {seq1_file}, {seq2_file}")
    print(f"Output prefix: {args.output}")
    
    # Perform alignment and get mappings
    A_to_B, B_to_A, seq1_id, seq2_id = perform_alignment(seq1_file, seq2_file, args.output)
    
    # Query specific site if requested
    if args.site:
        query_site_mapping(A_to_B, B_to_A, seq1_id, seq2_id, args.site)
    
    # Print sample mappings
    if not args.site:
        print(f"\nSample mappings (first 10 positions):")
        print(f"A to B mapping:")
        for i, (pos, mapping) in enumerate(sorted(A_to_B.items())):
            if i >= 10:
                break
            res_a, pos_b, res_b = mapping
            print(f"  A{pos}({res_a}) -> B{pos_b}({res_b})")
        
        print(f"B to A mapping:")
        for i, (pos, mapping) in enumerate(sorted(B_to_A.items())):
            if i >= 10:
                break
            res_b, pos_a, res_a = mapping
            print(f"  B{pos}({res_b}) -> A{pos_a}({res_a})")

if __name__ == "__main__":
    main()