#!/usr/bin/env python3
"""
PDB to FASTA Sequence Extractor
This script extracts protein sequences from PDB files and saves them in FASTA format.
"""

import os
import sys
from Bio import PDB
import argparse

def extract_sequences_from_pdb(pdb_file):
    """
    Extract all protein sequences from a PDB file
    Returns: dict {chain_id: sequence}
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)

        # Get the first model (most PDB files have only one model)
        model = structure[0]

        # Standard amino acid three-letter to one-letter mapping
        aa_dict = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        chain_sequences = {}

        for chain in model:
            chain_id = chain.id
            sequence = ""
            residue_count = 0

            for residue in chain:
                # Only consider standard amino acids (not water, ions, ligands, etc.)
                if residue.id[0] == ' ':  # Standard amino acid
                    res_name = residue.resname

                    if res_name in aa_dict:
                        one_letter = aa_dict[res_name]
                        sequence += one_letter
                        residue_count += 1
                    else:
                        # Non-standard amino acid - warn but continue
                        print(f"  Warning: Non-standard residue {res_name} at position {residue.id[1]} in chain {chain_id}")

            # Only add chains with protein sequences (at least 1 amino acid)
            if sequence:
                chain_sequences[chain_id] = sequence
                print(f"  Chain {chain_id}: {len(sequence)} residues")
            else:
                print(f"  Chain {chain_id}: No protein sequence found (skipped)")

        if not chain_sequences:
            print(f"  Warning: No protein sequences found in {pdb_file}")

        return chain_sequences

    except FileNotFoundError:
        print(f"Error: File {pdb_file} not found")
        return {}
    except Exception as e:
        print(f"Error reading PDB file {pdb_file}: {e}")
        return {}

def write_fasta_file(chain_sequences, pdb_file, output_file):
    """
    Write sequences to FASTA file
    """
    try:
        pdb_basename = os.path.splitext(os.path.basename(pdb_file))[0]

        with open(output_file, 'w') as f:
            for chain_id, sequence in chain_sequences.items():
                # Create FASTA header: >PDB_BASENAME_CHAIN_ID
                header = f">{pdb_basename}_{chain_id}"
                f.write(f"{header}\n")

                # Write sequence with line breaks every 80 characters
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

        print(f"  FASTA file written: {output_file}")
        return True

    except Exception as e:
        print(f"Error writing FASTA file {output_file}: {e}")
        return False

def process_pdb_file(pdb_file):
    """
    Process a single PDB file
    """
    print(f"\nProcessing: {pdb_file}")

    # Extract sequences
    chain_sequences = extract_sequences_from_pdb(pdb_file)

    if not chain_sequences:
        print(f"  No sequences extracted from {pdb_file}")
        return {}

    # Generate output filename
    pdb_basename = os.path.splitext(pdb_file)[0]
    output_file = f"{pdb_basename}.fasta"

    # Write FASTA file
    success = write_fasta_file(chain_sequences, pdb_file, output_file)

    if success:
        # Print chain:sequence dictionary
        print(f"  Chain:Sequence dictionary:")
        for chain_id, sequence in chain_sequences.items():
            # Print first 50 characters of sequence for readability
            seq_preview = sequence[:50] + "..." if len(sequence) > 50 else sequence
            print(f"    '{chain_id}': '{seq_preview}' (length: {len(sequence)})")

    return chain_sequences

def main():
    parser = argparse.ArgumentParser(description='Extract protein sequences from PDB files and convert to FASTA format')
    parser.add_argument('-i', '--input', nargs='+', required=True, metavar='PDB_FILE',
                       help='One or more PDB files to process')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Check if input files exist
    missing_files = []
    for pdb_file in args.input:
        if not os.path.exists(pdb_file):
            missing_files.append(pdb_file)

    if missing_files:
        print("Error: The following files were not found:")
        for file in missing_files:
            print(f"  - {file}")
        sys.exit(1)

    print("PDB to FASTA Converter")
    print("=" * 50)
    print(f"Processing {len(args.input)} PDB file(s)...")

    # Process all PDB files
    all_results = {}
    successful_files = 0

    for pdb_file in args.input:
        chain_sequences = process_pdb_file(pdb_file)
        if chain_sequences:
            all_results[pdb_file] = chain_sequences
            successful_files += 1

    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Total files processed: {len(args.input)}")
    print(f"Successful extractions: {successful_files}")
    print(f"Failed extractions: {len(args.input) - successful_files}")

    if all_results:
        print(f"\nOverall Results:")
        for pdb_file, chain_sequences in all_results.items():
            pdb_basename = os.path.basename(pdb_file)
            print(f"\n{pdb_basename}:")
            print(f"  Chains found: {list(chain_sequences.keys())}")
            print(f"  Total chains: {len(chain_sequences)}")

            if args.verbose:
                print(f"  Detailed sequences:")
                for chain_id, sequence in chain_sequences.items():
                    print(f"    Chain {chain_id}: {sequence}")

    # Print all chain:sequence dictionaries
    i#!/usr/bin/env python3
"""
PDB to FASTA Sequence Extractor
This script extracts protein sequences from PDB files and saves them in FASTA format.
"""

import os
import sys
from Bio import PDB
import argparse

def extract_sequences_from_pdb(pdb_file):
    """
    Extract all protein sequences from a PDB file
    Returns: dict {chain_id: sequence}
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)

        # Get the first model (most PDB files have only one model)
        model = structure[0]

        # Standard amino acid three-letter to one-letter mapping
        aa_dict = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        chain_sequences = {}

        for chain in model:
            chain_id = chain.id
            sequence = ""
            residue_count = 0

            for residue in chain:
                # Only consider standard amino acids (not water, ions, ligands, etc.)
                if residue.id[0] == ' ':  # Standard amino acid
                    res_name = residue.resname

                    if res_name in aa_dict:
                        one_letter = aa_dict[res_name]
                        sequence += one_letter
                        residue_count += 1
                    else:
                        # Non-standard amino acid - warn but continue
                        print(f"  Warning: Non-standard residue {res_name} at position {residue.id[1]} in chain {chain_id}")

            # Only add chains with protein sequences (at least 1 amino acid)
            if sequence:
                chain_sequences[chain_id] = sequence
                print(f"  Chain {chain_id}: {len(sequence)} residues")
            else:
                print(f"  Chain {chain_id}: No protein sequence found (skipped)")

        if not chain_sequences:
            print(f"  Warning: No protein sequences found in {pdb_file}")

        return chain_sequences

    except FileNotFoundError:
        print(f"Error: File {pdb_file} not found")
        return {}
    except Exception as e:
        print(f"Error reading PDB file {pdb_file}: {e}")
        return {}

def write_fasta_file(chain_sequences, pdb_file, output_file):
    """
    Write sequences to FASTA file
    """
    try:
        pdb_basename = os.path.splitext(os.path.basename(pdb_file))[0]

        with open(output_file, 'w') as f:
            for chain_id, sequence in chain_sequences.items():
                # Create FASTA header: >PDB_BASENAME_CHAIN_ID
                header = f">{pdb_basename}_{chain_id}"
                f.write(f"{header}\n")

                # Write sequence with line breaks every 80 characters
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

        print(f"  FASTA file written: {output_file}")
        return True

    except Exception as e:
        print(f"Error writing FASTA file {output_file}: {e}")
        return False

def process_pdb_file(pdb_file):
    """
    Process a single PDB file
    """
    print(f"\nProcessing: {pdb_file}")

    # Extract sequences
    chain_sequences = extract_sequences_from_pdb(pdb_file)

    if not chain_sequences:
        print(f"  No sequences extracted from {pdb_file}")
        return {}

    # Generate output filename
    pdb_basename = os.path.splitext(pdb_file)[0]
    output_file = f"{pdb_basename}.fasta"

    # Write FASTA file
    success = write_fasta_file(chain_sequences, pdb_file, output_file)

    if success:
        # Print chain:sequence dictionary
        print(f"  Chain:Sequence dictionary:")
        for chain_id, sequence in chain_sequences.items():
            # Print first 50 characters of sequence for readability
            seq_preview = sequence[:50] + "..." if len(sequence) > 50 else sequence
            print(f"    '{chain_id}': '{seq_preview}' (length: {len(sequence)})")

    return chain_sequences

def main():
    parser = argparse.ArgumentParser(description='Extract protein sequences from PDB files and convert to FASTA format')
    parser.add_argument('-i', '--input', nargs='+', required=True, metavar='PDB_FILE',
                       help='One or more PDB files to process')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Check if input files exist
    missing_files = []
    for pdb_file in args.input:
        if not os.path.exists(pdb_file):
            missing_files.append(pdb_file)

    if missing_files:
        print("Error: The following files were not found:")
        for file in missing_files:
            print(f"  - {file}")
        sys.exit(1)

    print("PDB to FASTA Converter")
    print("=" * 50)
    print(f"Processing {len(args.input)} PDB file(s)...")

    # Process all PDB files
    all_results = {}
    successful_files = 0

    for pdb_file in args.input:
        chain_sequences = process_pdb_file(pdb_file)
        if chain_sequences:
            all_results[pdb_file] = chain_sequences
            successful_files += 1

    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Total files processed: {len(args.input)}")
    print(f"Successful extractions: {successful_files}")
    print(f"Failed extractions: {len(args.input) - successful_files}")

    if all_results:
        print(f"\nOverall Results:")
        for pdb_file, chain_sequences in all_results.items():
            pdb_basename = os.path.basename(pdb_file)
            print(f"\n{pdb_basename}:")
            print(f"  Chains found: {list(chain_sequences.keys())}")
            print(f"  Total chains: {len(chain_sequences)}")

            if args.verbose:
                print(f"  Detailed sequences:")
                for chain_id, sequence in chain_sequences.items():
                    print(f"    Chain {chain_id}: {sequence}")

    # Print all chain:sequence dictionaries
    # if all_results:
    #     print(f"\nAll Chain:Sequence Dictionaries:")
    #     print("-" * 50)
    #     for pdb_file, chain_sequences in all_results.items():
    #         pdb_basename = os.path.basename(pdb_file)
    #         print(f"\n{pdb_basename}:")
    #         for chain_id, sequence in chain_sequences.items():
    #             if args.verbose:
    #                 print(f"  '{chain_id}': '{sequence}'")
    #             else:
    #                 seq_preview = sequence[:50] + "..." if len(sequence) > 50 else sequence
    #                 print(f"  '{chain_id}': '{seq_preview}' (length: {len(sequence)})")


if __name__ == "__main__":
    main()
