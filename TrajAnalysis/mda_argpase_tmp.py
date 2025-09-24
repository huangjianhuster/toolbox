# Author: Jian Huang
# coding=utf-8
# @Time :   12/23/22 7:59 AM
# @Software :   PyCharm
# @contact  :   jianhuang@umass.edu

import MDAnalysis
import numpy as np
import sys
import math
import argparse

### Global variables


### Functions




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count water number within a specific cone region')
    parser.add_argument('-s', '--structure', help='structure of the MD system (pdb format)', required=True)
    parser.add_argument('-t', '--trajectory', help='trajectory of the MD system (xtc format)')
    parser.add_argument('-i', '--resid', help='residue id of a selected potassium', required=True)
    parser.add_argument('-b', '--begin', help='the beginning frame to do counting', type=int, default=0)
    parser.add_argument('-e', '--end', help='the ending frame to do counting', type=int, default=999999)
    parser.add_argument('-k', '--skip', help='skip between frames for counting', type=int, default=1)

    args = parser.parse_args()
    structure_file = args.structure
    target_trajectory = args.trajectory
    initial_frame = args.begin
    final_frame = args.end
    skip = args.skip

    # load system
    if target_trajectory:
        target_system = MDAnalysis.Universe(structure_file, target_trajectory)
        # print("Length of the input trajectory", len(target_system.trajectory))
    else:
        target_system = MDAnalysis.Universe(structure_file)

   ### Your calculations below

