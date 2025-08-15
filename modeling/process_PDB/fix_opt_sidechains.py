import os
import sys
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.optimizers import ConjugateGradients

# Global variables
INPUT_PDB_FILE = sys.argv[1]
OUTPUT_PDB_FILE = sys.argv[2] if len(sys.argv) == 3 else INPUT_PDB_FILE.split('.')[0] + '_fixed.pdb'


### Phase 1: fix missing side chains
env = Environ()  # Create a new Modeller environment

env.io.atom_files_directory = ['.']  # Search for input PDB in the current directory
env.io.hetatm = False                # Do not read HETATM records from input PDB
env.io.water = False                 # Do not read WATER records from input PDB

env.libs.topology.read(file='$(LIB)/top_heav.lib') # Standard heavy atom topology library
env.libs.parameters.read(file='$(LIB)/par.lib')   # Standard parameter library

mdl = complete_pdb(env, filename=INPUT_PDB_FILE, special_patches=None, transfer_res_num=False, model_segment=None)

mdl.write(file=OUTPUT_PDB_FILE)

# Step 2: Select sidechains (to rebuild/refine)
atmsel = Selection(mdl).only_sidechain()

# Step 3: Create stereochemical restraints
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

# Optional: Randomize sidechains slightly to avoid local minima
atmsel.randomize_xyz(deviation=4.0)
# Step 4: Run energy minimization using conjugate gradients
cg = ConjugateGradients()
cg.optimize(atmsel, max_iterations=200)

# Step 5: Save the final model
mdl.write(file=OUTPUT_PDB_FILE)
