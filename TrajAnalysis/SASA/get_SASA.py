import os
import numpy as np
import mdtraj as md
from multiprocessing import Pool, cpu_count

def get_available_threads():
    """Get the number of available threads in the current environment."""
    total_threads = os.cpu_count()
    used_threads = psutil.cpu_count(logical=False)  # Physical cores
    available_threads = total_threads - used_threads
    return available_threads

def compute_sasa(traj_chunk):
    return md.shrake_rupley(traj_chunk, mode='residue') * 100


# Testing
if __name__ == "__main__":
    pdb = "./apo/trjs/eq_ProCen.pdb"
    xtc = "./apo/trjs/298Kr1_TMalign_dt1ns.xtc"

    # get available cores/threadings for parallel computing
    available_cores = get_available_threads()
    if n_cores > available_cores:
        n_cores = available_cores

    u = mda.Universe(pdb, xtc)
    traj = md.load(xtc, top=pdb)

    # sasa usually takes a long time; only select the part of the system you are interested
    # in mdtraj, the algorithm is aware of all atoms defined in the selection (traj)
    selection = traj.topology.select('protein or resname POPC')
    traj_with_lipid = traj.atom_slice(selection)

    chunk_size = len(traj_with_lipid) // n_cores
    chunks = [traj_with_lipid[i:i+chunk_size] for i in range(0, len(traj), chunk_size)]

    with Pool(n_cores) as pool:
        results = pool.map(compute_sasa, chunks)

    sasa_with_lipid = np.vstack(results)
    # shape: (n_frames, n_residues)

    # print(sasa_with_lipid)
    # np.savetxt("test.txt", sasa_with_lipid, fmt="%.2f", delimiter=',')
