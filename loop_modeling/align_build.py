# Author: Jian Huang
# coding=utf-8
# @Time :   9/16/22 5:54 PM
# @File :   align_build.py.py
# @Software :   PyCharm
# @contact  :   jianhuang@umass.edu

import os
from modeller import *
from modeller.automodel import *


# restrain all atom symmetry
class mymodel(automodel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A', 'B', 'C', 'D'], renumber_residues=[150, 150, 150, 150])
        self.chains[0].name = 'A'
        self.chains[1].name = 'B'
        self.chains[2].name = 'C'
        self.chains[3].name = 'D'

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


env = environ()
# directories for input atom files
env.io.atom_files_directory = './:../atom_files'
# Be sure to use 'mymodel' rather than 'automodel' here!
# knowns should be the code for the structure in the given alin file [the prefix of file name without .pdb]
# sequence should be the code for the sequence in the given alin file
a = mymodel(env, alnfile='8FC9_alignto8FC7.ali',
              knowns='8FC9_ABCD', sequence='8FC7_ABCD',
              assess_methods=(assess.DOPE,
                              assess.GA341))

a.starting_model = 1                # index of the first model
a.ending_model  = 5                # index of the last model
                                  # (determines how many models to calculate)
a.make()                           # do homology modeling
