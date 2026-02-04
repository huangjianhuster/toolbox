This README aims at:
1. explain scripts inside the `vmd` folder
2. exemplify the usages of those scripts using files from `./examples`


## `vmdrc`
This is my personalized vmdrc (initialization file for vmd), which defines the default settings for view, lights, colors, materials and other visualization settings.

To use the file

```bash
# save a copy of your current .vmdrc
mv ~/.vmdrc ~/.vmdrc_backup
cp ./vmdrc ~/.vmdrc
```

Below are some customized features:

1. Override the default `red2`, 'red3', 'blue2', 'blue3', 'green2', 'green3', 'orange2' and 'orange3' through re-defining their rgb values.
2. Set the default representation style as `AOChalky`+ `NewCartoon`, and color by `Name` with the selection of `protein`.
3. Fix window locations on the screen: check `menu main move` or `menu graphics move` (adjust those numbers according to your screen)
4. New materials: `Shiva`, `Fading` and `FadingSurf`
5. Utils functions:
    - `hide_rep`: hide all representations in the current "top" molecule
    - `show_rep`: show all representations in the current "top" molecule
    - `toggle_rep`: turn on or off a specific representation
    - `clean`: remove all representations in the current "top" molecules
    - `rep`: generate default representations (through predefined selections) for the current "top" molecule
    - `set_allcartoon`: modify parameters (thickness, resolution, aspect_ratio and spline style) for all Newcartoon representation (on the "top" molecule)
    - `vmd_draw_arrow`: drawing arrows


## `draw_COM.tcl`
Draw the center of mass of a selection with frame tracking.
- The default selection could be modified inside the script: `set selection "[your own selection]"`

```bash
# 0. load your molecule into VMD

# 1. load the script
source draw_COM.tcl
```

## `draw_interactions.tcl`
Draw pairwise residue interactions with cylinders in VMD graphic drawing.

Usage:
```bash
# 0. load your structure into VMD: vmd ./examples/8FC8.pdb
# 0. prepare your interaction data file with the following format:

# ChainID1  ResID1  ResName1  ChainID2  ResID2  ResName2
B           579  ILE       D           469  CYS     
B           355  VAL       B           503  HIS     
B           527  LYS       B           574  ASN     
B           350  TRP       B           292  HIS

or with interaction strength

# ChainID1  ResID1  ResName1  ChainID2  ResID2  ResName2  Strength
B           579  ILE       D           469  CYS       0.6592
B           355  VAL       B           503  HIS       0.1224
B           527  LYS       B           574  ASN       0.4368
B           350  TRP       B           292  HIS       0.5825



# 1. load tcl script
source draw_interaction.tcl

# 2. utils functions
# 2.1 print out basic structure and residue indices information
list_structure_info
# 2.2 read interaction data from the data file prepared (from above)
set interactions [read_data "./examples/interactions_wStrength.dat"]
# 2.3 draw interactions
draw_interactions $interactions

```

## `draw_label.tcl`
Draw labels of selected residues with frame tracking. This can be useful when visualizing contacting residues around a specific selection.

Usage:
```bash
# 0. load your structure and trajectory into VMD

# 1. load the script
source draw_label.tcl
```

You may need to modify the selection in the script (Line: 12) for your specifc use:

```tcl
set reslabel [atomselect top "resid 450 to 750 and name CA"]
```


## `alignment.tcl`
Perform alignment between two structures and save out selected atoms/atogroups.

This can be useful when you want to:
1. align two protein structures and save out the coordinates for further modeling or visualization
2. align two protein structures (one of which has ligands; for example align A to B, and A has ligands) and save out ligand-bound structure of B.

Usage:
```bash
# 0. load two structures into VMD, for example
# vmd -m A.pdb B.pdb 
# in this case, A.pdb will be the "move" molecule; B.pdb will be the "target" molecule in the following context

# 1. load the script
source alignment.tcl

# 2. performan alignment
# 2.1 for cases you want the same selection for alignment
align_and_save_ligand -align "protein and backbone" -ligand "all" -output "test.pdb"
# 2.2 for case you want different selections for the "move" and "target" molecules
align_and_save_ligand -align_move "chain A and name CA" \
                      -align_target "chain B and name CA" \
                      -ligand "all" \
                      -output "test.pdb"

```



## `analyze_structure.tcl`
Basic and simple statical analysis of system composition.

Usage:
```bash
# 0. load structure into VMD

# 1. load the script
source analyze_structure.tcl

# 2. run analysis
analyze_structure
list_chains
```

## `time.tcl`
Draw the time information of frames with frame tracking.

Usage:
```bash
source time.tcl
```

You may need to change the time step definition according to your trajectory saving frequency: `set timeperframe 0.05` (Line 26).
