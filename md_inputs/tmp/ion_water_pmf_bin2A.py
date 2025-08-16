import os
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
import math
import sys

# gloable variables
# pdb = "../step6.6.center.pdb"
# xtc = "../s1_TMalign.xtc"   # based on the TM aligned trajectory
pdb = sys.argv[1]
xtc = sys.argv[2]
output_suffix = sys.argv[3]

# define positions: center at gate (z=0)
u = mda.Universe(pdb,xtc)
gate_residue_394 = u.select_atoms("protein and resid 394 and name CA")
gate_residue_398 = u.select_atoms("protein and resid 398 and name CA")
gate_residue_394to398 = u.select_atoms("protein and resid 394:398 and backbone")
selectivity_filter = u.select_atoms("protein and resid 358:361 and name CA")
upper_memb = u.select_atoms("resname POPC and name P and prop z > 95")
lower_memb = u.select_atoms("resname POPC and name P and prop z < 95")
upper_memb_pos = upper_memb.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
lower_memb_pos = lower_memb.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
SF_pos = selectivity_filter.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
print("upper_memb_pos (A): ", upper_memb_pos)
print("lower_memb_pos (A)", lower_memb_pos)
print("selectivity filter (A)", SF_pos)

# calculation params
box_zsize = 140
x0, y0, z0 = gate_residue_394to398.center_of_mass()
bin_max = box_zsize - z0
bin_min = -1*z0
print(f"bin_min: {bin_min} Ang")
print(f"bin_max: {bin_max} Ang")

upper_bound = box_zsize     # end at largest box size
lower_bound = z0 + bin_min  # start from 0
K_z = []    # store K+ z positions
O_z = []    # store water O positions
frame_number = 0    # get whole frames
total_selected_K = 0    # counts total
total_selected_O = 0
cutoff = 4      # cutoff radius for the central cylinder

# calculation -- counts
for ts in u.trajectory:
    # Better to input a ns-timestep trajectory!
    # current_time = ts.time/1000.0 # in ns
    current_time = ts.time / 1000.0
    print("current_time", current_time)

    frame_number += 1
    # select potassium near pore region
    selection_K = "( resname POT )  and prop x >= %s and "  \
                "prop x <= %s and prop y >= %s and prop y <= %s " % (str(round(x0-cutoff)), str(round(x0+cutoff)),
                                                                        str(round(y0-cutoff)), str(round(y0+cutoff)))
        # select water molecules along the center pore region
    selection_O = "byres  (resname TIP3 and name OH2 and prop x >= %s and "  \
                "prop x <= %s and prop y >= %s and prop y <= %s )" % (str(round(x0-cutoff)), str(round(x0+cutoff)),
                                                                        str(round(y0-cutoff)), str(round(y0+cutoff)))
    # print(selection)
    potassium = u.select_atoms(selection_K)
    water_oxygen = u.select_atoms(selection_O)
    for K_atom_position in potassium.positions:
        x = float(K_atom_position[0])
        y = float(K_atom_position[1])
        z = float(K_atom_position[2])
        # print(z)
        # distance between K ion and the COM of gate residues
        if ((x-x0)*(x-x0)+(y-y0)*(y-y0) <= cutoff*cutoff) and (z >= lower_bound) and (z <= upper_bound):
            K_z.append(z - z0)
            total_selected_K += 1
    for O_atom_position in water_oxygen.positions:
        x = float(O_atom_position[0])
        y = float(O_atom_position[1])
        z = float(O_atom_position[2])
        # print(z)
        # distance between K ion and the COM of gate residues
        if (x-x0)*(x-x0)+(y-y0)*(y-y0) <= cutoff*cutoff:
            O_z.append(z - z0)
            # O_z.append(z)
            # print(z)
            total_selected_O += 1

# output for ion
def get_pmf(z_data, bin_num, bin_min, bin_max, saved_name=None):
    global frame_number, cutoff
    bins=bin_num
    hist, bin_edges = np.histogram(z_data, bin_num, (bin_min, bin_max))
    print("hist, bin_edges:", hist, bin_edges)
    sum_for_normalize = sum(hist)
    zaxis = np.array(bin_edges[:-1])
    grid = zaxis[1] - zaxis[0]
    counts = np.array(hist)
    
    kB = 1.3806488E-23  # j/K
    An = 6.02214179E23
    T = float(300)
    R = 1.987  # (cal/mol.degree)
    DG = np.zeros(bins)
    for x in range(bins):
        if counts[x] == 0:
            DG[x] = 99
        else:
            DG[x] = -0.001 * R * T * (math.log(((float(counts[x])) / frame_number) / (3.1415926 * cutoff * cutoff)))
            # 0.001 change cal to kcal
            # probability: counts/total frames/area; shouldn't it be Volume?
    out = np.vstack((zaxis, DG)).T
    if saved_name:
        np.savetxt(saved_name, out, header="Z-position(Ang)    PMF(kcal/mol)")
    return out


def get_density(z_data, bin_num, bin_min, bin_max, saved_name=None):
    global frame_number, cutoff
    bins=bin_num
    hist, bin_edges = np.histogram(z_data, bin_num, (bin_min, bin_max))
    print("hist, bin_edges:", hist, bin_edges)
    sum_for_normalize = sum(hist)
    zaxis = np.array(bin_edges[:-1])
    grid = zaxis[1] - zaxis[0]
    counts = np.array(hist)
    density = counts / (frame_number * 3.1415926 * cutoff * cutoff * grid) 
    out = np.vstack((zaxis, density)).T
    if saved_name:
        np.savetxt(saved_name, out, header="Z-position(Ang)    density")
    return out

# bin_min, bin_max
ion_pmf = get_pmf(K_z, bin_num=70, bin_min=bin_min, bin_max=bin_max, saved_name=f"{output_suffix}_dG_ion_bin2A.dat")
ion_density = get_density(K_z, bin_num=70, bin_min=bin_min, bin_max=bin_max, saved_name=f"{output_suffix}_ion_density_b2A_cf4.dat")
water_density = get_density(O_z, bin_num=70, bin_min=bin_min, bin_max=bin_max, saved_name=f"{output_suffix}_water_density_b2A_cf4.dat")
