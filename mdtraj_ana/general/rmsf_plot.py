# Author: Jian Huang
# Date: 2023/07/06

# TODO:
# plot RMSF data from PDB and rmsf calculation result;
# The PDB is used for secondard structure annotation

# Usage: 
# python rmsf_plot.py --pdb [PDB] --xvgfiles [1.xvg 2.xvg ...] --sep

# Dependencies
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import argparse
import warnings
import subprocess
warnings.filterwarnings("ignore")

# plt.style.use('/home2/jianhuang/projects/packages/toolbox/plot/mystyle.mplstyle')

# Functions
def run_stride(pdb_file):
    """Run STRIDE and return the output."""
    process = subprocess.run(["stride", pdb_file], capture_output=True, text=True)
    return process.stdout

def parse_stride_output(stride_output):
    """Parse the STRIDE output to get secondary structure ranges."""
    helices = []
    sheets = []
    loops = []

    for line in stride_output.splitlines():
        # Skip non-data lines
        if not line.startswith("LOC"):
            continue

        # Split STRIDE output line into columns
        columns = line.split()
        ss_type = columns[1]
        chainid = columns[4]
        ss_range = (int(columns[3]), int(columns[6]))

        if ss_type in ['AlphaHelix', '310Helix']:
            helices.append((chainid, ss_range))
        elif ss_type in ['Strand', ]:
            sheets.append((chainid, ss_range))
        else:
            loops.append((chainid, ss_range))
    return helices, sheets, loops

def get_data(xvg_files_list, average=True):
    ys = []
    labels = []
    for file in xvg_files_list:
        x, y = np.loadtxt(file,comments=["@","#"],unpack=True)
        # y *= 10
        ys.append(y)
        label = os.path.basename(file).split('.')[0]
        labels.append(label)
    if average:
        y_mean = np.mean(ys, axis=0)
        y_std = np.std(ys, axis=0)
        return x, y_mean, y_std
    else:
        return x, ys, labels

def main():
    # user defined variables
    parser = argparse.ArgumentParser(description="plot RMSF figure",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb", help="pdb file", required=True)
    parser.add_argument("--xvgfiles", help="gmx-xvg files for rmsf calculation", required=True, nargs="+")
    parser.add_argument("--sep", help="plot separate rmsf plots; otherwise plot average rmsf", action='store_true')
    args = parser.parse_args()

    # get secondary structures
    pdb = args.pdb
    resid, ss = get_ss(pdb)
    # Simplify SS results
    stride_output = run_stride(pdb_file=pdb)
    helix_ranges, sheet_ranges, loop_ranges = parse_stride_output(stride_output)
    
    # plot function
    fig, ax = plt.subplots(1,1,figsize=(10,6))

    # rmsf plot with average and errorbars
    rmsf_xvg_files = args.xvgfiles
    if not args.sep:
        x, y_mean, y_std = get_data(rmsf_xvg_files, average=True) 
        ax.errorbar(x, y_mean, yerr=y_std, elinewidth=1, linewidth=1.5, capsize=2, label='Average')
        y_max = max(y_mean) + (max(y_mean) - min(y_mean))*0.2
    else: # plot separately
        x, ys, labels = get_data(rmsf_xvg_files, average=False)
        y_max = 0
        for y,label in zip(ys, labels):
            ax.plot(x, y, "-", label=label)
            y_max_tmp = max(y) + (max(y) - min(y))*0.2
            if y_max < y_max_tmp:
                y_max = y_max_tmp
    
    # add secondary structure information
    for chainid,(start,end) in helix_ranges:
        ax.fill_betweenx([0,y_max], start, end, color='lightpink', alpha=0.4)
    for chainid,(start,end) in sheet_ranges:
        ax.fill_betweenx([0,y_max], start, end, color='lightgreen', alpha=0.4)
    
    # ax.grid()
    ax.set_xlim(min(x), max(x))
    ax.set_ylim([0, y_max])
    ax.set_ylabel("RMSF")
    ax.set_xlabel("Res ID")
    ax.legend()

    plt.show()

if __name__ == "__main__":
    
    main()
