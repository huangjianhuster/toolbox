import os
import sys
from pprint import pprint
from modeller import *
from modeller.automodel import *

# Input variable (PDB path)
INPUT_PDB_FILE = sys.argv[1]

# Modeller
# log.verbose()    # Request verbose output (useful for debugging)
env = Environ()  # Create a new Modeller environment

env.io.atom_files_directory = ['.']  # Search for input PDB in the current directory
env.io.hetatm = False                # Do not read HETATM records from input PDB
env.io.water = False                 # Do not read WATER records from input PDB

# Read the backbone PDB file
mdl = Model(env, file=INPUT_PDB_FILE)

aln = Alignment(env)
input_code = INPUT_PDB_FILE.split('.')[0]
aln.append_model(mdl, align_codes=input_code, atom_files=INPUT_PDB_FILE)

structure = aln[0]
seq = {}
for chain in structure.chains:
    seq[chain.name] = ''
    for residue in chain.residues:
        seq[chain.name] += residue.code
pprint(seq)
