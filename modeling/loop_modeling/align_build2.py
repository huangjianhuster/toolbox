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
class mymodel(LoopModel):
    # you can renumber the chain ID and residue number ID
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A', 'B', 'C', 'D'], renumber_residues=[1, 1, 1, 1])
        self.chains[0].name = 'A'
        self.chains[1].name = 'B'
        self.chains[2].name = 'C'
        self.chains[3].name = 'D'

    # symmetry constraints
    # def special_restraints(self, aln):
    #     # Constrain the A and B chains to be identical (but only restrain
    #     # the C-alpha atoms, to reduce the number of interatomic distances
    #     # that need to be calculated):
    #     s1 = selection(self.residue_range('459:A', '467:A'), self.residue_range('997:A', '1007:A'))
    #     s2 = selection(self.residue_range('459:B', '467:B'), self.residue_range('997:B', '1007:B'))
    #     s3 = selection(self.residue_range('459:C', '467:C'), self.residue_range('997:C', '1007:C'))
    #     s4 = selection(self.residue_range('459:D', '467:D'), self.residue_range('997:D', '1007:D'))
    #     self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
    #     self.restraints.symmetry.append(symmetry(s2, s3, 1.0))
    #     self.restraints.symmetry.append(symmetry(s3, s4, 1.0))

    # def user_after_single_model(self):
    #     # Report on symmetry violations greater than 1A after building
    #     # each model:
    #     self.restraints.symmetry.report(1.0)

    # Selected atoms to restraint the modeling region to be the missing loop
    # ref: https://salilab.org/modeller/manual/node23.html#SECTION:model-segment
    def select_loop_atoms(self):
        # Select residue 30 ~ 40 in chain A
        return Selection(
                         self.residue_range('343:A', '348:A'),
                         self.residue_range('459:A', '467:A'),
                         self.residue_range('989:A', '1000:A'),
                         self.residue_range('343:B', '348:B'),
                         self.residue_range('459:B', '467:B'),
                         self.residue_range('989:B', '1000:B'),
                         self.residue_range('343:C', '348:C'),
                         self.residue_range('459:C', '467:C'),
                         self.residue_range('989:C', '1000:C'),
                         self.residue_range('343:D', '348:D'),
                         self.residue_range('459:D', '467:D'),
                         self.residue_range('989:D', '1000:D'),
                         )


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
a.ending_model  = 1                # index of the last model
                                  # (determines how many models to calculate)

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 3           # Last loop model
a.loop.refine = 'fast'
a.make()                           # do homology modeling
