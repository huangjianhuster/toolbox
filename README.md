# Introduction
a generic toolbox for computational biophysicists


## Visualization
`dir: VMD`:
- an example `.vmdrc` for default settings
- `morph.tcl`: morphing between two structures of a protein
- `rotate_z.tcl`: rotate a protein (with coordinates change) in vmd
- `draw_COM.tcl`: draw COMs of a atomgroup selection in vmd
- `time.tcl`: draw time in vmd
- `draw_label.tcl`: draw labels in vmd
- `draw_interactions`: draw interactions in vmd
- `secondary_structure_cache.tcl`: cache secondary structure calculation in vmd
- `draw_principalaxes.tcl`: draw principal axes for a rigid body (domains)

`PyMol`:
- `pymolrc`: my own pymolrc settings

## MD trajectory analysis
`dir: TrajAnalysis`:
- `mda_argpase_tmp.py`: a template for argument control
- `mda_threading_example.py`: a template for threading your analysis
- `channel`: for ion channel analysis
- `general`: for general proteins 
- `SecondaryStructureAnalysis`: calculating secondary structure over a trajectory

## Modeling
- `dir: loop_modeling`: model missing loops using 'Modeller'
- `dir: alignment`: generate sequence alignment from fasta files or pdb files
- `dir: process_PDB`: 
    - extract sequence from PDB files: `get_seq_from_pdb.py`. Dependencies: Biopython
    - extract frames from xtc: `extract_frames_from_xtc.py`. Dependencies: MDAnalysis
    - model side chains (assuming all backbone atoms are present): `fix_opt_sidechains.py`. Dependencies: Modeller


## Scientific plot
`plot`: some of my own matplotlib practices and examples. Check the [plot practices](./plot/README.md)

## Bash functions
`bash/utils`: some useful bash functions for file I/O and parallelization.
