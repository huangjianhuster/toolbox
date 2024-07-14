import numpy as np
import matplotlib.pyplot as plt

# Data path <-- already removed PBC
pot1_positions_rmPBC = np.loadtxt("./analysis/s1_POT1_pos_rmPBC.dat", skiprows=1)
pot2_positions_rmPBC = np.loadtxt("./analysis/s1_POT2_pos_rmPBC.dat", skiprows=1)
pot3_positions_rmPBC = np.loadtxt("./analysis/s1_POT3_pos_rmPBC.dat", skiprows=1)
pot4_positions_rmPBC = np.loadtxt("./analysis/s1_POT4_pos_rmPBC.dat", skiprows=1)
pot5_positions_rmPBC = np.loadtxt("./analysis/s1_POT5_pos_rmPBC.dat", skiprows=1)

# Plot
fig, ax = plt.subplots(1,1,figsize=(16,8))
ax.plot(pot1_positions_rmPBC[:, 0], pot1_positions_rmPBC[:, 1])
ax.plot(pot2_positions_rmPBC[:, 0], pot2_positions_rmPBC[:, 1])
ax.plot(pot3_positions_rmPBC[:, 0], pot3_positions_rmPBC[:, 1])
ax.plot(pot4_positions_rmPBC[:, 0], pot4_positions_rmPBC[:, 1])
ax.plot(pot5_positions_rmPBC[:, 0], pot5_positions_rmPBC[:, 1])

# ax.set_ylim([-15,15])
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel("Simulation Time (ns)", fontsize=16)
ax.set_ylabel("Z-position (ns)", fontsize=16)
ax.set_xlim([-10,510])

# ax.axvline(16)
ax.axhspan(-5, 30, alpha=0.3, color='gray')
# ax.legend()
# plt.savefig("./500mV_SFhaHalfKcal_S61kcal/ion_permeation_location.svg")
# plt.show() 
 
 
# Plot -- zoom in
fig2, ax2 = plt.subplots(1,1,figsize=(16,6))
x_range = pot1_positions_rmPBC[:, 0]
ax2.plot(x_range[0:250], pot1_positions_rmPBC[0:250, 1], label="ID-1716", alpha=1, lw=3, zorder=2)
ax2.plot(x_range[:1200], pot2_positions_rmPBC[:1200, 1], label="ID-1530", alpha=1, lw=3, zorder=3)
ax2.plot(x_range[1400:2550], pot3_positions_rmPBC[1400:2550, 1], label="ID-1264", alpha=1, lw=3, zorder=4)
ax2.plot(x_range[:4600], pot4_positions_rmPBC[:4600, 1], label="ID-1447", alpha=0.7, lw=3, zorder=1)
ax2.plot(x_range[5000:7400], pot5_positions_rmPBC[5000:7400, 1], label="ID-1531", alpha=1, lw=3, zorder=5)

ax2.set_ylim([-10,35])
ax2.set_xlim([-10,510])
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.set_xlabel("Simulation Time (ns)", fontsize=16)
ax2.set_ylabel(r"Z-position ($\AA$)", fontsize=16)

# ax.axvline(16)
ax2.axhline(30.77, c='k', ls='dashed', label="upper membrane")
ax2.axhline(-7.98, c='k', ls='dashed',label="lower membrane")
ax2.grid()
# ax.legend(fontsize=16)
plt.savefig("./analysis/permeation.svg")
plt.show()
