import MDAnalysis as mda
from MDAnalysis.coordinates import DCD
import numpy as np
import argparse
import string
import subprocess

def fix_chain_ids_with_convpdb(infile, outfile):
    cmd = ["convpdb.pl", "-chainfromseg", "-out", "generic", infile]
    with open(outfile, "w") as fout:
        subprocess.run(cmd, stdout=fout, check=True)

def reassign_chain_ids(atomgroup, outname):
    # Assign sequential chain IDs (A, B, C, ...)
    chain_ids = string.ascii_uppercase  # 'A' to 'Z'
    segs = atomgroup.segments

    if len(segs) > len(chain_ids):
        raise ValueError("More segments than available chain IDs (A-Z). Extend mapping manually.")

    for seg, cid in zip(segs, chain_ids):
        seg.segid = 'PRO' + cid  # update segid to chain ID

    # Write to PDB
    atomgroup.write(outname)
    print(f"Saved with chain IDs reassigned → {outname}")

def main(psf_file, xtc_file, selection="backbone", output_prefix="random_frames",
         n_frames=3, seed=42, frames=None):
    u = mda.Universe(psf_file, xtc_file)
    selected = u.select_atoms(selection)

    total_frames = len(u.trajectory)
    if total_frames == 0:
        raise RuntimeError("Trajectory has zero frames.")

    if frames:  # user-specified frames
        frame_indices = [int(f) for f in frames]
        print(f"Using user-defined frames: {frame_indices}")
    else:  # random selection
        rng = np.random.default_rng(seed)
        frame_indices = rng.choice(total_frames, size=n_frames, replace=False)
        print(f"Randomly selected frames: {frame_indices}")
    # Write out selected frames
    for i, frame_idx in enumerate(frame_indices, start=1):
        u.trajectory[frame_idx]  # move to the frame
        raw_name = f"{output_prefix}_{frame_idx}_raw.pdb"
        fixed_name = f"{output_prefix}_{frame_idx}.pdb"
        reassign_chain_ids(selected, raw_name)
        fix_chain_ids_with_convpdb(raw_name, fixed_name)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract random backbone frames from trajectory.")
    parser.add_argument("psf", help="Input PSF file")
    parser.add_argument("xtc", help="Input XTC file")
    parser.add_argument("--prefix", default="random_frames", help="Output file prefix")
    parser.add_argument("--nframes", type=int, default=3, help="Number of random frames")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--frames", nargs="+", help="User-defined frame indices (space-separated)")
    parser.add_argument("--selection", default="backbone", help="mdanalysis selection syntax")

    args = parser.parse_args()
    main(args.psf, args.xtc, args.selection, args.prefix, args.nframes, args.seed, args.frames)
