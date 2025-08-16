import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

# matplotlib.rcParams['font.family'] = 'Arial'

# global variables
ion_pmfs = ["./analysis/s1__dG_ion_bin2A.dat",
            "./analysis/s2__dG_ion_bin2A.dat",
            "./analysis/s3__dG_ion_bin2A.dat"]

ion_densities = ["./analysis/s1__ion_density_b2A_cf4.dat",
                "./analysis/s2__ion_density_b2A_cf4.dat",
                "./analysis/s3__ion_density_b2A_cf4.dat"]

water_densities = ["./analysis/s1__water_density_b2A_cf4.dat",
                "./analysis/s2__water_density_b2A_cf4.dat",
                "./analysis/s3__water_density_b2A_cf4.dat",
                ]

# water coordination number
water_coor_datafile = "/home2/jianhuang/projects/hcn/open_refined/500mV_SFhaHalfKcal_S61kcal/coor_water/water_coordination.dat"
off_center_datafile = "/home2/jianhuang/projects/hcn/open_refined/500mV_SFhaHalfKcal_S61kcal/coor_water/off_center.dat"
water_coor_data = np.loadtxt(water_coor_datafile, skiprows=1)
off_center_data = np.loadtxt(off_center_datafile, skiprows=1)



upper_memb_pos, lower_memb_pos, SF_pos = 30.77, -7.98, (21.7671, 28.732)

# get ion pmf average and std
def get_ave_std(file_list):
    data_list = []
    for i in file_list:
        tmp = np.loadtxt(i, skiprows=1)
        data_list.append(tmp)
    data_array = np.concatenate(data_list, axis=1)
    data_array = data_array[:,[0,1,3,5]]
    x = data_array[:,0]
    ave = data_array[:,1:].mean(axis=1)
    std = data_array[:,1:].std(axis=1)
    return x, ave, std


# pre-processing for water coordination and offset of the pore axis data
def get_ave_std2(data, bin_edges=np.arange(-90, 60, 2)):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # Bin the positions (first column)
    positions = data[:, 0]
    hist, _ = np.histogram(positions, bins=bin_edges)

    # Calculate the mean and standard deviation of the second column for each bin
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    averages = []
    std_devs = []
    for i in range(len(bin_edges) - 1):
        bin_mask = (positions >= bin_edges[i]) & (positions < bin_edges[i + 1])
        values_in_bin = data[bin_mask, 1]
        if len(values_in_bin) > 0:
            averages.append(np.mean(values_in_bin))
            std_devs.append(np.std(values_in_bin))
        else:
            averages.append(np.nan)  # If no values in the bin, assign NaN
            std_devs.append(np.nan)

    return bin_centers, averages, std_devs


x_pmf, ion_pmf_ave, ion_pmf_std = get_ave_std(ion_pmfs)
x_ion_density, ion_density_ave, ion_density_std = get_ave_std(ion_densities)
x_water_density, water_density_ave, water_density_std = get_ave_std(water_densities)

bin_centers_water_coor, averages_water_coor, std_devs_water_coor = get_ave_std2(water_coor_data)
bin_centers_off_center, averages_off_center, std_devs_off_center = get_ave_std2(off_center_data)



fig, ax = plt.subplots(1,4,figsize=(16,6))
# pmf_ion = np.loadtxt(ion_pmf, skiprows=1)
# ax[0].plot(pmf_ion.T[1]-pmf_ion.T[1][-1], pmf_ion.T[0], "go-", label="K+ PMF")
ax[0].errorbar(ion_pmf_ave-ion_pmf_ave[-1], x_pmf, xerr=ion_pmf_std, capsize=2, elinewidth=1, fmt="g-", lw=2)
ax[0].set_xlabel("PMF(kcal/mol)", fontsize=16, )
ax[0].set_xlim([-2,5])
ax[0].set_title("Ion PMF", fontname="Arial")

# density_dat = np.loadtxt(ion_density, skiprows=1)
# ax[1].plot(density_dat.T[1], density_dat.T[0], "bo-", label="K+ density")
# ax[1].errorbar(ion_density_ave, x_ion_density, xerr=ion_density_std, capsize=2, elinewidth=1, fmt="b-", lw=2)
# ax[1].set_xlabel("Density", fontsize=14)
# ax[1].set_title("Ion Density")

# density_dat_water = np.loadtxt(water_density, skiprows=1)
# ax[2].plot(density_dat_water.T[1] / density_dat_water.T[1][-1], density_dat_water.T[0], "co-", label="water density")
# print("water_density_std", water_density_std)
ax[1].errorbar(water_density_ave/water_density_ave[-1], x_water_density, xerr=water_density_std/water_density_ave[-1], capsize=2, elinewidth=1, fmt="c-", lw=2)
ax[1].set_xlabel("Density", fontsize=14, fontname="Arial")
ax[1].set_title("Water Density", fontname="Arial")

ax[2].errorbar(averages_water_coor, bin_centers_water_coor, xerr=std_devs_water_coor, fmt='-', color='tomato', capsize=2, elinewidth=1,lw=2)
ax[2].set_xlabel('Water coordination number', fontsize=14, fontname="Arial")
# ax[2].set_ylabel(r'Z-position($\AA$)', fontsize=14)
ax[2].set_xlim([1,10])
ax[2].set_xticks([1, 2, 4, 6, 8, 10])
ax[2].set_title("Water Coordination Number", fontname="Arial")

ax[3].errorbar(averages_off_center[40:60], bin_centers_off_center[40:60], xerr=std_devs_off_center[40:60], fmt='b-', capsize=2, elinewidth=1,lw=2)
ax[3].errorbar(averages_off_center[:41], bin_centers_off_center[:41], xerr=std_devs_off_center[:41], fmt='b-', capsize=2, elinewidth=1,lw=2, alpha=0.2)
ax[3].errorbar(averages_off_center[59:], bin_centers_off_center[59:], xerr=std_devs_off_center[59:], fmt='b-', capsize=2, elinewidth=1,lw=2, alpha=0.2)
ax[3].set_xlabel('Ion offset from the pore axis', fontsize=14, fontname="Arial")
# ax[3].set_ylabel(r'Z-position($\AA$)', fontsize=14)
ax[3].set_xlim([0,14])
ax[3].set_xticks([0, 4, 8, 12, 14])
ax[3].set_title("Distance to the pore axis", fontname="Arial")


for i in ax.flatten():
    i.set_ylabel("Z-axis", fontsize=14, fontname="Arial")
    i.tick_params(axis='both', which='major', labelsize=14)
    i.grid()
    # annotations
    i.set_ylim([-75, 45])
    i.axhline(upper_memb_pos, color='orange', linestyle="--")
    i.axhline(lower_memb_pos, color='orange', linestyle="--")
    i.axhline(0, color='red', linestyle="--")
    i.axhspan(ymin=SF_pos[0], ymax=SF_pos[1], color='gray', alpha=0.5, linestyle="--")








plt.tight_layout()
plt.savefig("./analysis/pmf_ave2.svg")
plt.show()
