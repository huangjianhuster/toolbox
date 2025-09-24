import os
import math
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from multiprocessing import Pool
from functools import partial
import sys

# Usage:
# python pore_wat_cal.py [PDB] [XTC] [OUT]

### Selection Functions ###
def water_selection(x0, y0, cutoff):
    """
    rough region defined by x = x0 +/- cutoff; y = y0 +/- cutoff
    """
    selection_syntax = "byres  (resname TIP3 and name OH2 and prop x >= %s and "  \
                "prop x <= %s and prop y >= %s and prop y <= %s )" % (str(round(x0-cutoff)), str(round(x0+cutoff)),
                                                                        str(round(y0-cutoff)), str(round(y0+cutoff)))
    return selection_syntax

# select potassium ions
def pot_selection(x0, y0, cutoff):
    """
    rough region defined by x = x0 +/- cutoff; y = y0 +/- cutoff
    """
    selection_syntax = " (resname POT) and prop x >= %s and "  \
                "prop x <= %s and prop y >= %s and prop y <= %s" % (str(round(x0-cutoff)), str(round(x0+cutoff)),
                                                                        str(round(y0-cutoff)), str(round(y0+cutoff)))
    return selection_syntax


### Theading Calculation ###
# Base function for threading
def get_counts(frame_index, box_xy_center, cutoff, z_shift, ags):
    global z_lower_bound, z_upper_bound
    # ags is supposed to be a list of atomgroups
    ag1 = ags[0]
    ag1.universe.trajectory[frame_index]

    x0, y0 = box_xy_center

    result = []
    for ag in ags:
        counts = []
        # total_counts = 0
        for pos in ag.positions:
            x, y, z = pos
            if ((x-x0)*(x-x0)+(y-y0)*(y-y0) <= cutoff*cutoff) and (z >= z_lower_bound) and (z <= z_upper_bound):
                counts.append(z - z_shift)
                # total_counts += 1
        result.append(counts)

    return result

# Multiprocessing threading
def get_counts_threading(u, box_xy_center, radius, z_shift, atomgroups, n_jobs=6):

    # multiple groups could happen here
    # run_per_frame = partial(get_counts, box_xy_center=box_xy_center, cutoff=radius, z_shift=z_shift, ags=[group1, ]) 
    run_per_frame = partial(get_counts, box_xy_center=box_xy_center, cutoff=radius, z_shift=z_shift, ags=atomgroups) 
    
    frame_values = np.arange(u.trajectory.n_frames)

    with Pool(n_jobs) as worker_pool:
        result = worker_pool.map(run_per_frame, frame_values)
    return result

def get_box_xy_center(pdb_file):
    u_tmp = mda.Universe(pdb_file)
    x_dim, y_dim, *_ = u_tmp.dimensions
    return x_dim/2, y_dim/2

def get_box_z_dim(pdb_file):
    u_tmp = mda.Universe(pdb_file)
    z_dim = u_tmp.dimensions[2]
    return z_dim


### Data Export Function ###
def get_density(z_data, bin_num, bin_min, bin_max, frame_number, cutoff, saved_name=None):
    hist, bin_edges = np.histogram(z_data, bin_num, (bin_min, bin_max))
    # print("hist, bin_edges:", hist, bin_edges)
    zaxis = np.array(bin_edges[:-1])
    grid = zaxis[1] - zaxis[0]
    counts = np.array(hist)
    density = counts / (frame_number * 3.1415926 * cutoff * cutoff * grid) 
    out = np.vstack((zaxis, density)).T
    if saved_name:
        np.savetxt(saved_name, out, fmt='%.3f',  header="Z-position(Ang)    density")
        return None
    else:
        return out

def get_pmf(z_data, bin_num, bin_min, bin_max, frame_number, cutoff, saved_name=None):
    hist, bin_edges = np.histogram(z_data, bin_num, (bin_min, bin_max))
    # print("hist, bin_edges:", hist, bin_edges)
    sum_for_normalize = sum(hist)
    zaxis = np.array(bin_edges[:-1])
    grid = zaxis[1] - zaxis[0]
    counts = np.array(hist)
    
    kB = 1.3806488E-23  # j/K
    An = 6.02214179E23
    T = float(300)
    R = 1.987  # (cal/mol.degree)
    DG = np.zeros(bin_num)
    for x in range(bin_num):
        if counts[x] == 0:
            DG[x] = 99
        else:
            DG[x] = -0.001 * R * T * (math.log(((float(counts[x])) / frame_number) / (3.1415926 * cutoff * cutoff)))
            # 0.001 change cal to kcal
            # probability: counts/total frames/area; shouldn't it be Volume?
    out = np.vstack((zaxis, DG)).T
    if saved_name:
        np.savetxt(saved_name, out, fmt='%.3f', header="Z-position(Ang)    PMF(kcal/mol)")
        return None
    else:
        return out

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


def main(pdb_file, xtc_file, radius, z_shift, n_jobs=8):
    u = mda.Universe(pdb_file, xtc_file)
    x_center, y_center = get_box_xy_center(pdb_file)
    print(f"box xy dimension center: {x_center}\t{y_center}")
    
    group1 = u.select_atoms(water_selection(x0=x_center, y0=y_center, cutoff=radius), updating=True)
    atomgroups = [group1, ]
    result = get_counts_threading(u, box_xy_center=(x_center, y_center), radius=radius, z_shift=z_shift, atomgroups=atomgroups, n_jobs=n_jobs)

    # print(result[0])
    return result

if __name__ == '__main__':
    # path to PDB and trajectories
    pdb = sys.argv[1]
    xtc = sys.argv[2]
    out = sys.argv[3]
  
    # Global variables
    # by default, calculate the whole box
    z_lower_bound = 0
    z_upper_bound = get_box_z_dim(pdb_file=pdb)
    # default radius = 4 angstrom to select the cylindral region
    cutoff=4
    # by default, no z_shift
    z_shift=0
    bins = int((z_upper_bound-z_lower_bound)*1.5)

    data = main(pdb_file=pdb, xtc_file=xtc, radius=cutoff, z_shift=z_shift)
    # data will be a list of list, each sublist is atom_z_coordinates of a specific frame
    # data[frame_num][0]: coordinates of atomgroup1 of the frame_num
    # data[frame_num][1]: coordiantes of atomgroup2 of the frame_num
    # ...

    # 1 angstrom as the bin size
    # here, the code only calculate water density, so only 1 kind of atomgroup
    # <-- it is possible to include ions; add atomgroups in to the varible in the main function
    group_num = 1
    group_name = ['wat', ]

    frame_number = len(data)
    for i,name in zip(range(group_num), group_name):
        unpack_tmp=[]
        for frame in data:
            unpack_tmp.extend([ag_coords for idx, ag_coords in enumerate(frame) if idx==i])
        unpack_arr = np.concatenate(unpack_tmp)
        get_density(z_data=unpack_arr, bin_num=bins, bin_min=z_lower_bound, bin_max=z_upper_bound,\
                            frame_number=frame_number, cutoff=cutoff, saved_name=f"{out}_{name}_density.dat")
        get_pmf(z_data=unpack_arr, bin_num=bins, bin_min=z_lower_bound, bin_max=z_upper_bound,\
                            frame_number=frame_number, cutoff=cutoff, saved_name=f"{out}_{name}_pmf.dat")

