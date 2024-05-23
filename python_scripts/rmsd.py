# Author: Jian Huang
# E-mail: jianhuang@umass.edu
# TODO: calculate RMSD 

import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import sys

# Load input
PSF = sys.argv[1]
TRJ = sys.argv[2]
u = MDAnalysis.Universe(PSF, TRJ)
ref = MDAnalysis.Universe(PSF,TRJ)     # reference closed AdK (1AKE) (with the default ref_frame=0)

# Align for superimpose; 
# of note, gromacs backbone doesnot have Oxygen
align_selection = "protein and backbone and not name O"


R = rms.RMSD(u, ref,
           select=align_selection,             # superimpose on whole backbone of the whole protein
           groupselections=["backbone",        # CORE
                            ])                 # NMP
R.run()


