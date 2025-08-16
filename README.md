# toolbox
a toolbox for Computational Biophysics

## Loop modeling

dir: `loop_modeling`

Goal: Fixing missing loops in a PDB file.

## VMD scripts

dir: `vmd_scripts`

Goal: tcl scripts used in VMD visualization. 

Usage: (open VMD)

```text
source xxx.tcl
```

- `vmdrc`: vmd initiation file. Replace your '~/.vmdrc' with this file.
- `time.tcl`: add simulation time to the OpenGL window of VMD.
- `morph.tcl`: morphing between two PDB structures.
- `draw_interactions.tcl`: draw interacting residue pairs.
- `default_rep.tcl`: quickly generate representations for your system.
- `default_rep_memb.tcl`; quickly generate representations for membrane protein system.

### 
