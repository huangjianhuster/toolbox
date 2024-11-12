# Author: Jian Huang
# coding=utf-8
# @Time :   6/18/23 
# @File :   align_input.py
# @Software :   PyCharm
# @contact  :   jianhuang@umass.edu


from modeller import *
import sys

# generate original alin file for a pdb
pdb_file = sys.argv[1]
pdb_code = pdb_file.split('.')[0]
aln_file = pdb_file.split('.')[0] + ".aln"

e = Environ()
aln = Alignment(e)
m = Model(e, file=pdb_file)
aln.append_model(m, align_codes=pdb_code, atom_files=pdb_file)
aln.write(file=aln_file)
