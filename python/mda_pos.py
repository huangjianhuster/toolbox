import MDAnalysis as mda
import sys

# gloable variables
pdb = sys.argv[1]

# define positions: center at gate (z=0)
u = mda.Universe(pdb)
gate_residue_394 = u.select_atoms("protein and resid 394 and name CA")
gate_residue_398 = u.select_atoms("protein and resid 398 and name CA")
gate_residue_394to398 = u.select_atoms("protein and resid 394:398 and backbone")
selectivity_filter_358 = u.select_atoms("protein and resid 358 and name CA")
selectivity_filter = u.select_atoms("protein and resid 358:361 and name CA")
selectivity_filter_361 = u.select_atoms("protein and resid 361 and name CA")
upper_memb = u.select_atoms("resname POPC and name P and prop z > 95")
lower_memb = u.select_atoms("resname POPC and name P and prop z < 95")

upper_memb_pos = upper_memb.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
lower_memb_pos = lower_memb.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
SF_pos = selectivity_filter.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
SF_lower_358_pos = selectivity_filter_358.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]
SF_upper_361_pos = selectivity_filter_361.center_of_mass()[2] - gate_residue_394to398.center_of_mass()[2]

print("upper_memb_pos (A): ", upper_memb_pos)
print("lower_memb_pos (A): ", lower_memb_pos)
print("selectivity filter (A): ", SF_pos)
print("SF 358 (A): ", SF_lower_358_pos)
print("SF 361 (A): ", SF_upper_361_pos)

