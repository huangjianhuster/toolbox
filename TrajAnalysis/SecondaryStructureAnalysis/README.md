# Secondary structure calculation

There are generally two popular methods for secondary structure assignment: **DSSP** and **STRIDE**.

Check the following webpages for their publications:

```
# DSSP
Kabsch W, Sander C (1983). “Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features”. Biopolymers 22 (12): 2577-637. doi:10.1002/bip.360221211

# STDIE
Frishman D, Argos P. [Knowledge-Based Protein Secondary Structure Assignment](https://webclu.bio.wzw.tum.de/stride/stride.pdf) Proteins: Structure, Function, and Genetics 23:566-579 (1995)
```

## Use DSSP 

DSSP is implemented in both the `mdtraj` and `MDAnalysis`.

1. DSSP in `MDAnalysis`: follow the [tutorial page](https://docs.mdanalysis.org/dev/documentation_pages/analysis/dssp.html)
2. DSSP in `mdtraj`:

```bash
python mdtraj_dssp.py [PDB] [XTC] [OUT_PREFIX]
```

Modify the script to only calculate secondary structures for a specific selection as explained in the `mdtraj_dssp.py` header.

## Use STRIDE

STRIDE is the default internal secondary structure calculation algorithm implemented in VMD. The provided `vmd_stride.tcl` can help.

Modify the arguments written in the `vmd_stride.tcl` to run your calculations, such as PDB, xtc, selections etc.

When runing the script, we can mutethe visualization of VMD to speed up the calculation and reduce memory usage:

```bash
vmd -dispdev text -e vmd_stride.tcl
```

An example output would be like the following:

```text
Frame   B:205   B:206   B:207   B:208   B:209   B:210
0       H       H       H       C       C       C
1       H       H       H       C       C       C
2       H       H       H       H       T       T
3       H       H       H       H       C       C
4       T       T       T       T       T       T
5       H       H       T       T       T       C
6       H       T       T       T       T       C
...
...
45      H       H       H       C       C       C
46      H       H       H       H       C       C
47      H       H       H       H       C       C
48      H       H       T       T       T       C
49      H       H       H       H       C       C
50      H       H       H       H       C       C
51      H       H       H       H       C       C
```

The secondary structure annotations (such as 'H', 'T', 'C') follow the STRIDE convention:

```text
H	    Alpha helix
G	    3-10 helix
I	    PI-helix
E	    Extended conformation
B or	b   Isolated bridge
T	    Turn
C	    Coil (none of the above)
```

