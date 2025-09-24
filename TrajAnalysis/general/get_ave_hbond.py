import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import numpy as np
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


def get_ave_hbond_tyr386(u):
    """
    u: MDanalysis. Universe (load trjectory)
    """
    # selection of the Y386
    tyr_selection = "protein and resname TYR and resid 223 686 1149 1612 and not backbone"
    
    hbonds= HBA(universe=u, donors_sel=tyr_selection, acceptors_sel=tyr_selection, d_a_cutoff=3.5)
    # default cutoff angle: 150
    hbonds.run()
    
    ave_hbond = hbonds.results['hbonds'].shape[0] / 4 / u.trajectory.n_frames
    return ave_hbond


# load trajectory
# paths to trajectories
pdb = "/home2/jianhuang/projects/hcn/open_refined/500mV_S61kcal/step6.6_equilibration.center.pdb"
tpr = "/home2/jianhuang/pikes_home/work/hcn/open_refined/production/step6.6_equilibration.tpr"

# S6 unrestrained; only 100 ns simulation
traj_S6_unrestrained_1 = "/home2/jianhuang/projects/hcn/open_refined/500mV_SF_0.5kcal/s1_500mV_SF_half.xtc"
traj_S6_unrestrained_1_part2 = "/home2/jianhuang/projects/hcn/open_refined/500mV_SF_0.5kcal/s1_100ns-500ns.part0002.xtc"
traj_S6_unrestrained_1_part3 = "/home2/jianhuang/projects/hcn/open_refined/500mV_SF_0.5kcal/s1_100ns-500ns.part0003.xtc"

u = mda.Universe(tpr, traj_S6_unrestrained_1, traj_S6_unrestrained_1_part2, traj_S6_unrestrained_1_part3)
ave_hbond = get_ave_hbond_tyr386(u)
print(f"Average hydrogen bond of Y386: {ave_hbond}")
