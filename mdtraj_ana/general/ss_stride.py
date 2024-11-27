import subprocess
import sys

# Dependencies: STRIED
# install: https://webclu.bio.wzw.tum.de/stride/install.html

# Usage:
# python ss_stride.py [PDB] 

def run_stride(pdb_file):
    """Run STRIDE and return the output."""
    process = subprocess.run(["stride", pdb_file], capture_output=True, text=True)
    return process.stdout

def parse_stride_output(stride_output):
    """Parse the STRIDE output to get secondary structure ranges."""
    helices = []
    sheets = []
    loops = []

    for line in stride_output.splitlines():
        # Skip non-data lines
        if not line.startswith("LOC"):
            continue

        # Split STRIDE output line into columns
        columns = line.split()
        ss_type = columns[1]
        chainid = columns[4]
        ss_range = (int(columns[3]), int(columns[6]))

        if ss_type in ['AlphaHelix', '310Helix']:
            helices.append((chainid, ss_range))
        elif ss_type in ['Strand', ]:
            sheets.append((chainid, ss_range))
        else:
            loops.append((chainid, ss_range))
    return helices, sheets, loops


def main(pdb_file):
    """Main function to calculate and parse secondary structures."""
    # Run STRIDE to get secondary structure data
    stride_output = run_stride(pdb_file)

    # Parse the STRIDE output
    helix_ranges, sheet_ranges, loop_ranges = parse_stride_output(stride_output)

    # Print the results
    for chainid,i in helix_ranges:
        if i[1]-i[0]>4:
            print(f"Helix: {chainid}-{i[0]}~{i[1]}")
    for chainid,i in sheet_ranges:
        if i[1]-i[0]>4:
            print(f"Sheet: {chainid}-{i[0]}~{i[1]}")


if __name__ == "__main__":
    pdb_file = sys.argv[1 ]# "eq_res_protein.pdb"  # Replace with your PDB file
    main(pdb_file)

