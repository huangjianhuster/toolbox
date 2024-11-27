import json
import matplotlib.pyplot as plt
import sys
import numpy as np

# load input json
filelist = sys.argv[1:]
data_list = []
for file in filelist:
    with open(file, "r") as f:
        data = json.load(f)
        data_list.append(data)


def extract_pore(data):
    # extract pore profile information
    pore_profile = data['pathwayProfile']
    pore_axis = np.array(pore_profile['s'])*10
    pore_radius_mean = np.array(pore_profile['radiusMean'])*10
    pore_radius_std = np.array(pore_profile['radiusSd'])*10

    return pore_axis, pore_radius_mean, pore_radius_std

# get pore data list
pore_data_list = []
for data in data_list:
    pore_data_list.append(extract_pore(data))

# plot
fig, ax = plt.subplots(1, 1, figsize=(6,10))
for pore_data in pore_data_list:
    pore_axis, pore_radius_mean, pore_radius_std = pore_data[0], pore_data[1], pore_data[2]
    ax.errorbar(pore_radius_mean, pore_axis, xerr=pore_radius_std, fmt="-", capsize=2, lw=2, elinewidth=1, errorevery=10)

ax.set_xlim([0,16])
ax.set_ylim([-75, 45])
ax.set_xlabel(r"Pore radius ($\AA$)", fontsize=14)
ax.set_ylabel(r"Z-axis ($\AA$)", fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid()
plt.tight_layout()
plt.savefig("pore-radius.svg")
plt.show()
