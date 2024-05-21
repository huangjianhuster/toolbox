# Author: Jian Huang
# Date: 2023/07/06

# TODO:
# plot rmsf

# Dependencies
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from Bio.PDB import PDBParser
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import warnings
warnings.filterwarnings("ignore")

# plt.style.use('/home2/jianhuang/projects/packages/toolbox/plot/mystyle.mplstyle')

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

def get_data(xvg_files_list, average=True):
    ys = []
    for file in xvg_files_list:
        x, y = np.loadtxt(file,comments=["@","#"],unpack=True)
        # y *= 10
        ys.append(y)
    if average:
        y_mean = np.mean(ys, axis=0)
        y_std = np.std(ys, axis=0)
        return x, y_mean, y_std
    else:
        return x, ys

def main():
    # user defined variables
    parser = argparse.ArgumentParser(description="plot RMSF figure",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb", help="pdb file", required=True)
    parser.add_argument("--xvgfiles", help="gmx-xvg files for rmsf calculation", required=True, nargs="+")
    parser.add_argument("--sep", help="plot separate rmsf plots; otherwise plot average rmsf", action='store_true')
    parser.add_argument("--labels", help="labels", default=None, nargs="+")
    args = parser.parse_args()

    # get secondary structures
    pdb = args.pdb
    resid, ss = get_ss(pdb)
    # Simplify SS results
    for i in range(len(ss)):
        # print(ss[i])
        if ss[i] == 'G' or ss[i] == 'I': # make pi helix and 310 helix to be "helix"
            ss[i] = 'H'
        elif ss[i] == '-' or ss[i] == 'T' or ss[i] == 'S' or ss[i]== 'B':   # make other types into "C" (coil)
            # print(ss[i])
            ss[i] = 'C'
    print("\n-------Secondary Structure---------")
    # print(ss)
    
    # plot function
    fig, ax = plt.subplots(1,1,figsize=(10,6))

    # rmsf plot with average and errorbars
    rmsf_xvg_files = args.xvgfiles
    if not args.sep:
        x, y_mean, y_std = get_data(rmsf_xvg_files, average=True) 
        ax.errorbar(x, y_mean, yerr=y_std, elinewidth=1, linewidth=1.5, capsize=4)
        y_max = max(y_mean) + (max(y_mean) - min(y_mean))*0.2
    else: # plot separately
        x, ys = get_data(rmsf_xvg_files, average=False)
        y_max = 0
        for y in ys:
            ax.plot(x, y, "-")
            y_max_tmp = max(y) + (max(y) - min(y))*0.2
            if y_max < y_max_tmp:
                y_max = y_max_tmp
    
    # ax.grid()
    ax.set_xlim(min(x), max(x))
    ax.set_ylim([0, y_max])
    ax.set_ylabel("RMSF")
    ax.set_xlabel("Res ID")
    
    # add secondary structure information
    first = 0
    ss_range = []
    # print(ss)
    for i in range(len(ss)-1):
        # first = i
        last = i
        if ss[i] == ss[i+1]:
            last = i+1
            # print(last, len(ss))
    
        else:
            # fill area
            # x_a = [x[first], x[last]]
            # print(x[first], x[last])
            if ss[first] == 'H':
                ax.fill_betweenx([0,y_max], x[first], x[last], color='lightpink', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Helix")
            elif ss[first] == 'E':
                ax.fill_betweenx([0,y_max], x[first], x[last], color='lightgreen', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Sheet")
            else:
                ax.fill_betweenx([0,y_max], x[first], x[last], color='white', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Coil")
            first = i+1
            
        if last == len(ss)-1:
            # print(x[first], x[last])
            if ss[first] == 'H':
                ax.fill_betweenx([0,y_max], x[first], x[last], color='lightpink', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Helix")
            elif ss[first] == 'E':
                ax.fill_betweenx([0,y_max], x[first], x[last], color='lightgreen', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Sheet")
            else:
                ax.fill_betweenx([0,y_max], x[first], x[last], color='white', alpha=0.25)
                print(f"{int(x[first])}-{int(x[last])}: Coil")
    plt.show()

if __name__ == "__main__":
    
    main()
