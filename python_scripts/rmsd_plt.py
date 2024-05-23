import numpy as np
import matplotlib.pyplot as plt
import sys

# Load datafile
datafiles = sys.argv[1:]
data_arrays = []
for datafile in datafiles:
    data = np.loadtxt(datafile, comments=['@', '#'])
    data[:, 1] *= 10    # gromcas nm --> ang
    data_arrays.append(data)


# plot
fig, ax = plt.subplots(figsize=(12,4))
for data,name in zip(data_arrays, datafiles):
    ax.plot(data[:, 0]/1000, data[:, 1], lw=2, label=f"{name}")

ax.set_xlabel("Simulation time (ns)", fontsize=14)
ax.set_yticks([0,1,2,3,4,5,6,7,8])
ax.set_ylabel(r"RMSD ($\AA$)", fontsize=14)
ax.set_ylim([0,8])
ax.tick_params(axis='both', which='major', labelsize=14)

plt.legend(frameon=False)
plt.grid()
plt.tight_layout()
plt.savefig("output-rmsd.svg")
plt.show()


