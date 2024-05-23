import numpy as np
import matplotlib.pyplot as plt

# global variables
# ion_pmf = "./ion_dG_bin2A.dat"
# ion_density = "./ion_density_b1A_cf4.dat"
# water_density = "./water_b2A_cf4_f100ns.dat"
ion_pmf = sys.argv[1]
ion_density = sys.argv[2]
water_density = sys.argv[3]

upper_memb_pos, lower_memb_pos, SF_pos = 29.216, -7.667, (20.527, 28.292)


fig, ax = plt.subplots(1,3,figsize=(12,6))
pmf_ion = np.loadtxt(ion_pmf, skiprows=1)
ax[0].plot(pmf_ion.T[1]-pmf_ion.T[1][-1], pmf_ion.T[0], "go-", label="K+ PMF")
ax[0].set_xlabel("PMF(kcal/mol)", fontsize=16)
ax[0].set_xlim([-1,10])
ax[0].set_title("Ion PMF")

density_dat = np.loadtxt(ion_density, skiprows=1)
ax[1].plot(density_dat.T[1], density_dat.T[0], "bo-", label="K+ density")
ax[1].set_xlabel("Density", fontsize=14)
ax[1].set_title("Ion Density")

density_dat_water = np.loadtxt(water_density, skiprows=1)
ax[2].plot(density_dat_water.T[1] / density_dat_water.T[1][-1], density_dat_water.T[0], "co-", label="water density")
ax[2].set_xlabel("Density", fontsize=14)
ax[2].set_title("Water Density")

for i in ax.flatten():
    i.set_ylabel("Z-axis", fontsize=14)
    i.tick_params(axis='both', which='major', labelsize=14)
    i.grid()
    # annotations
    i.axhline(upper_memb_pos, color='r', linestyle="--")
    i.axhline(lower_memb_pos, color='r', linestyle="--")
    i.axhline(0, color='orange', linestyle="--")
    i.axhspan(ymin=SF_pos[0], ymax=SF_pos[1], color='gray', alpha=0.5, linestyle="--")

plt.tight_layout()
plt.savefig("ion_water_pmf_f100ns.svg")
plt.show()
