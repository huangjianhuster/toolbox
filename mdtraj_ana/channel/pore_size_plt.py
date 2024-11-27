import matplotlib.pyplot as plt
import numpy as np
import os
import pickle 
import sys

# input data files including three cols: Z-axis, radius_mean, radius_std
# pickle_file = sys.argv[1]
allfiles = sys.argv[1:]

if os.path.exists("/home2/jianhuang/projects/packages/toolbox/plot/mystyle.mplstyle"):
    plt.style.use("/home2/jianhuang/projects/packages/toolbox/plot/mystyle.mplstyle")


def read_pickle_file(pickle_file):
    with open(pickle_file, "rb") as pickle_file:
        loaded_dict = pickle.load(pickle_file)
    for key in loaded_dict.keys():
        print(f"{key:<10}{len(loaded_dict[key]):>10}")
    return loaded_dict



def pore_dim_plt(files):
    plt.figure(figsize=(4,8))
    for avestd_file in files:
        avestd_data = np.loadtxt(avestd_file, comments="#")
    
        # plot
        z_coor = avestd_data[:, 0]
        radii_mean = avestd_data[:, 1]
        radii_std = avestd_data[:, 2]
        gate_z = 0
    
        plt.plot(radii_mean, z_coor-gate_z, '-')
        plt.fill_betweenx(z_coor-gate_z, radii_mean-radii_std, radii_mean+radii_std, alpha=0.2)
    
    # set labels
    plt.xlim(0, 16)
    plt.xticks([0,4,8,12,16])
    # plt.ylim(-35, 80)
    # plt.yticks([-35, -20, 0, 20, 40, 60, 80])
    plt.xlabel(r"Pore radius (${\AA}$)")
    plt.ylabel("Z-position (${\AA}$)")
    
    # annotation
    # plt.axhline(upper_memb, color='orange', linestyle="--")
    # plt.axhline(lower_memb, color='orange', linestyle="--")
    # plt.axhline(0, color='red', linestyle="--")
    # plt.axhspan(ymin=selectivity_filter_469, ymax=selectivity_filter_472, color='gray', alpha=0.3, linestyle="--")
    
    plt.tight_layout()
    # plt.savefig("pore_size_plt.svg")
    plt.show()


pore_dim_plt(allfiles)
