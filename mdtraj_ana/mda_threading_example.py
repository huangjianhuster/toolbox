import os
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from multiprocessing import Pool
from functools import partial
import sys

# Usage:
# python mda_threading_example.py [PDB] [XTC] [OUT]

# path to PDB and trajectories
pdb = sys.argv[1]
xtc = sys.argv[2]
out = sys.argv[3]

# Load mda universe
u = mda.Universe(pdb, xtc)

# This is an example to calculate the Z-directional COM distance between two atomgroups 
# as a function of simulation time

### Helper Functions ###
# SF+PH; selection function 1
def SF_selection(chainid):
    return f"protein and chainid {chainid} and resid 624-640 and backbone"


# S4: selection function 2
def S4_selection(chainid):
    return f"protein and chainid {chainid} and resid 546-566 and backbone"


# Calculate Z-directional COM distances
def get_zCOM(u, selection1, selection2):
    com1_z = u.select_atoms(selection1).center_of_mass()[-1]
    com2_z = u.select_atoms(selection2).center_of_mass()[-1]
    return com2_z-com1_z

### Theading Calculation ###
# Base function for threading
def get_zCOM_perframe(frame_index, ag1, ag2):
    ag1.universe.trajectory[frame_index]
    distance = ag1.center_of_mass()[-1] - ag2.center_of_mass()[-1]
    return ag1.universe.trajectory.time, distance

# Multiprocessing threading
def get_all_zCOM(u, chainid, n_jobs=6):
    group1 = u.select_atoms(SF_selection(chainid))
    group2 = u.select_atoms(S4_selection(chainid))

    run_per_frame = partial(get_zCOM_perframe, ag1=group1, ag2=group2)
    frame_values = np.arange(u.trajectory.n_frames)

    with Pool(n_jobs) as worker_pool:
        result = worker_pool.map(run_per_frame, frame_values)
    return np.asarray(result)


### Data Export Function ###
def save_data(*arr1, filename, header=None):
    """
    arr1: a list of 1D arrays; for manipulating all kinds of arrays
    """
    tmp = np.vstack(arr1).T
    if header:
        np.savetxt(filename, tmp, fmt='%.3f', delimiter='\t', header=header)
    else:
        np.savetxt(filename, tmp, fmt='%.3f', delimiter='\t')
    return None


chA = get_all_zCOM(u=u, chainid='A')
chB = get_all_zCOM(u=u, chainid='B')
chC = get_all_zCOM(u=u, chainid='C')
chD = get_all_zCOM(u=u, chainid='D')

save_data([chA[:,0], chA[:, 1], chB[:, 1], chC[:, 1], chD[:, 1]], filename=out,\
         header="#Time\tchain_A\tchain_B\tchain_C\tchain_D")

