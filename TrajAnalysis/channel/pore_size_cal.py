import MDAnalysis as mda
from MDAnalysis.analysis import hole2
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import warnings
import pickle
warnings.filterwarnings('ignore')

top = sys.argv[1]
xtc = sys.argv[2]
outdir = sys.argv[3]

if not os.path.exists("./tmp/"):
    os.mkdir("./tmp/")

if not os.path.exists(f"{outdir}"):
    os.mkdir(f"{outdir}")


# Load trajectory
u = mda.Universe(top, xtc)
total_frames = len(u.trajectory)

# Hole calculation
ha = hole2.HoleAnalysis(u, select='protein',
                        cpoint='center_of_geometry',
                        cvect=[0,0,1],
                        end_radius=15.0,
                        tmpdir="./tmp/",
                        prefix=f"./tmp/{xtc.split('.')[0]}",
                        # executable='~/programs/hole2/exe/hole',
                        )
ha.run()

# get radii and edges
radii, edges = ha.bin_radii(bins=100, range=None)
# get mean and std
radii_mean = np.array([i.mean() for i in radii])
radii_std = np.array([i.std() for i in radii])
midpoints = 0.5*(edges[1:]+edges[:-1])

# Save to a pickle file
with open(f"{outdir}/{xtc.split('.')[0]}.pkl", "wb") as pickle_file:
    pickle.dump(ha.gather(), pickle_file)

ha.delete_temporary_files()

def save_data(*arr1, filename, header=None):
    tmp = np.vstack(arr1).T
    if header:
        np.savetxt(filename, tmp, fmt='%.3f', delimiter='\t', header=header)
    else:
        np.savetxt(filename, tmp, fmt='%.3f', delimiter='\t')
    return None

save_data(midpoints, radii_mean, radii_std, filename=f"{outdir}//{xtc.split('.')[0]}.avestd.dat",\
         header="Z-axis\tradius_mean(Ang)\tradius_std")


