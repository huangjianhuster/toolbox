# Author: Jian Huang
# coding=utf-8
# @Time :   9/16/22 5:54 PM
# @File :   align_build.py.py
# @Software :   PyCharm
# @contact  :   jianhuang@umass.edu

import os
from modeller import *
from modeller.automodel import *
import sys

align_file = sys.argv[1]
knowns_code = sys.argv[2]
seq_code = sys.argv[3]

# restrain all atom symmetry
class mymodel(automodel):
    # you can renumber the chain ID and residue number ID
    """
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A', 'B', 'C', 'D'], renumber_residues=[1, 1, 1, 1])
        self.chains[0].name = 'A'
        self.chains[1].name = 'B'
        self.chains[2].name = 'C'
        self.chains[3].name = 'D'
    """

    # symmetry constraints
    """
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(self.chains['A'])
        s2 = selection(self.chains['B'])
        s3 = selection(self.chains['C'])
        s4 = selection(self.chains['D'])
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
        self.restraints.symmetry.append(symmetry(s2, s3, 1.0))
        self.restraints.symmetry.append(symmetry(s3, s4, 1.0))

    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)
    """

    # Selected atoms to restraint the modeling region to be the missing loop
    # ref: https://salilab.org/modeller/manual/node23.html#SECTION:model-segment
    def select_atoms(self):
        # Select residue 30 ~ 40 in chain A
        return Selection(self.residue_range('30:A', '40:A'))


env = environ()
# directories for input atom files
env.io.atom_files_directory = './:../atom_files'

# Be sure to use 'mymodel' rather than 'automodel' here!

# knowns should be the code for the structure in the given alin file [the prefix of file name without .pdb]
# sequence should be the code for the sequence in the given alin file

# get the code for knowns and sequence

a = mymodel(env, alnfile=align_file,
              knowns=knowns_code, sequence=seq_code,
              assess_methods=(assess.DOPE,
                              assess.GA341))

a.starting_model = 1                # index of the first model
a.ending_model  = 5                # index of the last model
                                  # (determines how many models to calculate)
a.make()                           # do homology modeling
