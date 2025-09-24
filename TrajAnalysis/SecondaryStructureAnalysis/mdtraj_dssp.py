import numpy as np
import mdtraj as md
import os
import sys
sys.path.append("/home2/jianhuang/projects/trpm3/scripts/")
from utils import *

# The following is an example to:
# calculate RESID 205 - 221 secondary structures for four monomers with chainID: B, F, J, N
# adjust the function accordingly for your own specific system

def get_ss_mdtraj(pdb, xtc, chains=['B', 'F', 'J', 'N']):
    traj = md.load(xtc, top=pdb)
    top = traj.topology
    # create a mapping from PDB-based chainID/residueID to mdtraj-based residueID
    # mdtraj by default numbers residues sequentially from 0. (0, 1, 2, ...); and ignore
    # chain information.
    pdb_map = {(res.chain.chain_id, res.resSeq): res for res in top.residues}
    ave_helicity = []
    ave_loop = []
    ave_sheet = []
    for chainid in chains:
        # chain B and resid 205 to 220
        sel_residues = [pdb_map[(chainid, r)] for r in range(205, 221)\
                        if (chainid, r) in pdb_map]

        # This is for debugging and checking whether we have selected correct residues in PDB
        # print out selected residues
        # for res in sel_residues:
        #     print(f"selected: ChainID-{res.chain.chain_id}:ResID-{res.name}{res.resSeq}")

        atom_indices = [atom.index for res in sel_residues for atom in res.atoms]

        # slice trajectory based on atom indices
        sub_traj = traj.atom_slice(atom_indices)
        dssp = md.compute_dssp(sub_traj, simplified=True)
        helicity_per_residue = np.mean(dssp == 'H', axis=0)

        loop_per_residue = np.mean(dssp == 'C', axis=0)
        sheet_per_residue = np.mean(dssp == 'E', axis=0)

        ave_helicity.append(helicity_per_residue)
        ave_loop.append(loop_per_residue)
        ave_sheet.append(sheet_per_residue)

    return np.array(ave_loop).mean(axis=0),\
           np.array(ave_helicity).mean(axis=0),\
           np.array(ave_sheet).mean(axis=0)

pdb = sys.argv[1] # "/home2/jianhuang/projects/trpm3/apo_18/rms/step4.1_equilibration.pdb"
xtc = sys.argv[2] # "/home2/jianhuang/projects/trpm3/apo_18/rms/r1_500ns_303K_cyto_dt1ns.xtc"
prefix = sys.argv[3]

results = get_ss_mdtraj(pdb=pdb, xtc=xtc)
save_ss(np.array(results).T, f"{prefix}.txt")
