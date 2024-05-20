# Author: Jian Huang
# Date: 2023/07/06

# TODO:
# plot rmsf

# USAGE:
# python rmsf_plot.py [PDB_file]

# multiple xvg files?
# average?
# annotate SS automatically

# Dependencies
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import sys
import argparse
# parse PDB
from Bio.PDB import PDBParser
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import sys

# Global variables
chainA_pdb = "/home2/jianhuang/projects/hydrophobic_dewetting/TRPV4/PIP2/RMSF/A.pdb"
chainE_pdb = "/home2/jianhuang/projects/hydrophobic_dewetting/TRPV4/PIP2/RMSF/E.pdb"

# Functions
class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

def get_ss(pdb_file):
    p = PDBParser()
    pdb = "pdb_file"
    structure = p.get_structure(' ', pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    resid = []
    ss = []
    for i,j in zip(list(dssp.keys()), list(model.get_residues())):
        # resid.append(dssp[i][0])
        ss.append(dssp[i][2])
        resid.append(j.full_id[-1][1])
    return resid, ss

def get_ave(xvg_files_list):
    ys = []
    for file in xvg_files_list:
        print(file)
        x,y = np.loadtxt(file,comments=["@","#"],unpack=True)
        y *= 10
        # reconver original numbering
        ys.append(y)
        x_new = [i+147 for i in range(1, 492+1)] + [i+147 for i in range(514, 639+1)]

        # mask array
        mc = ma.array(y)
        mc[492] = ma.masked
        # ax.plot(x_new, mc, label=os.path.basename(file), alpha=0.7)

    y_mean = np.mean(ys, axis=0)
    y_std = np.std(ys, axis=0)
    mc = ma.array(y_mean)
    mc[492] = ma.masked
    return x_new, mc, y_std

def rmsf_plt(xvg_files_list, color, ss):
    """
    ss: list of secondary structures
    """
    fig, ax = plt.subplots(figsize=(30, 5))
    ys = []
    for file in xvg_files_list:
        print(file)
        x,y = np.loadtxt(file,comments=["@","#"],unpack=True)
        y *= 10
        # reconver original numbering
        ys.append(y)
        x_new = [i+147 for i in range(1, 492+1)] + [i+147 for i in range(514, 639+1)]

        # mask array
        mc = ma.array(y)
        mc[492] = ma.masked
        # ax.plot(x_new, mc, label=os.path.basename(file), alpha=0.7)

    y_mean = np.mean(ys, axis=0)
    y_std = np.std(ys, axis=0)
    mc = ma.array(y_mean)
    mc[492] = ma.masked
    # ax.plot(x_new, mc, "k--", linewidth=3, label="ave")
    makers, caps, bars = ax.errorbar(x_new, mc, yerr=y_std, elinewidth=1, linewidth=2, capsize=2, color=color, ecolor=color, errorevery=2)
    # [bar.set_alpha(0.5) for bar in bars]
    
    ax.grid()
    ax.set_xlim(min(x_new), max(x_new))
    ax.set_xticks([min(x_new), 200, 300, 400, 500, 600, 639, 661, 786])
    ax.set_ylim([0, 7])
    ax.set_ylabel("RMSF(Ang)")
    ax.set_xlabel("Res ID")
    
    # annotations
    ax.plot([148, 404], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("ARD", (276, 0.4), weight='bold')
    ax.plot([320+147, 344+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S1", ((320+147+344+147)/2, 0.4), weight='bold')
    ax.plot([358+147, 387+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S2", ((358+147+387+147)/2, 0.4), weight='bold')
    ax.plot([400+147, 422+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S3", ((400+147+422+147)/2, 0.4), weight='bold')
    ax.plot([425+147, 447+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S4", ((425+147+447+147)/2, 0.4), weight='bold')
    ax.plot([455+147, 488+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S5", ((455+147+488+147)/2, 0.4), weight='bold')
    ax.plot([545+147, 574+147], [0.3, 0.3], 'k-', linewidth=3)
    ax.annotate("S6", ((545+147+574+147)/2, 0.4), weight='bold')
    
    # add secondary structure information
    first = 0
    ss_range = []
    for i in range(len(ss)-1):
        # first = i
        last = i
        if ss[i] == ss[i+1]:
            last = i+1
            # print(last, len(ss))

        else:
            # fill area
            # x_a = [x_new[first], x_new[last]]
            if ss[first] == 'H':
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightpink', alpha=0.25)
            elif ss[first] == 'E':
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightgreen', alpha=0.25)
            else:
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightyellow', alpha=0.25)
            first = i+1
            
        if last == len(ss)-1:
            if ss[first] == 'H':
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightpink', alpha=0.25)
            elif ss[first] == 'E':
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightgreen', alpha=0.25)
            else:
                ax.fill_betweenx([0,7], x_new[first], x_new[last], color='lightyellow', alpha=0.25)
    plt.show()
    return x_new, mc, y_std


def main():
    # user defined variables
    parser = argparse.ArgumentParser(description="plot RMSF figure",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--pdb", help="pdb file", required=True)
    parser.add_argument("-f", "--files", help="gmx-xvg files for rmsf calculation", required=True, nargs="+")
    parser.add_argument("-c", "--color", help="a plt-recognizable color", default='b', nargs="+")
    parser.add_argument("--sep", help="plot together or separate", action='store_true')
    parser.add_argument("-l", "--labels", help="labels", default=None, nargs="+")
    args = parser.parse_args()


    # get secondary structures
    pdb = args.pdb
    resid, ss = get_ss(pdb)
    
    for i in range(len(ss)):
        # print(ss[i])
        if ss[i] == 'G' or ss[i] == 'I': # make pi helix and 310 helix to be "helix"
            ss[i] = 'H'
        elif ss[i] == '-' or ss[i] == 'T' or ss[i] == 'S' or ss[i]== 'B':   # make other types into "C" (coil)
            # print(ss[i])
            ss[i] = 'C'

    # rmsf plot with average and errorbars
    list_of_rmsf_paths = args.files

    x, y, y_std = get_ave_std(list_of_rmsf_paths)




    for file in list_of_rmsf_paths:
        pass
    if args.sep:
        for i in list_of_rmsf_paths:
            x, y_ave, y_std = rmsf_plt(list_of_rmsf_paths[4*i: 4*(i+1)], color, ss)
    else:
        datalist = []
        for i , label in zip(range(num), args.labels):
            x, y_ave, y_std = get_ave(list_of_rmsf_paths[4*i: 4*(i+1)])
            datalist.append((x, y_ave, y_std, label))
        multiple_plt(datalist, ss, args.color)


if __name__ == "__main__":
    
    main()
