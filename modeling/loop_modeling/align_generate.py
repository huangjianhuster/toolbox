# Author: Jian Huang
# coding=utf-8
# @Time :   9/16/22 5:32 PM
# @File :   align_generate.py.py
# @Software :   PyCharm
# @contact  :   jianhuang@umass.edu


from modeller import *
import sys

# generate alin input file for a modeler
pdb_file = sys.argv[1]
seq_file = sys.argv[2]

# for example, this pdb structure has some missing loops
pdb_code = pdb_file.split('.')[0]
with open(seq_file, 'r') as file:
    first = file.readline().strip()
seq_code = first.split(';')[-1]

# making inputs
e = Environ()
aln = Alignment(e)
m = Model(e, file=pdb_file)

# add pdb
aln.append_model(m, align_codes=pdb_code, atom_files=pdb_file)

# add the new sequence you want to build [need to write by your own based your own needs!]
# example: in this alin file, your can add those missing loops
aln.append(file=seq_file, align_codes=seq_code)

aln.align2d(max_gap_length=50)
aln.write(file='seq_alignment.ali', alignment_format='PIR')
aln.write(file='seq_alignment.pap', alignment_format='PAP')
# remember to check the ali file to make the alignment is correct!!!
