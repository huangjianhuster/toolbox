# Python scripts for MD trajectory analysis and plot

## Template python files of using MDAnalysis
1. Mdanalysis + argument parser template: `mda_argpase_tmp.py`
2. Mdanalysis + Multiprocessing threading: `mda_threading_example.py`

## Generic analysis
dir: `./general`

1. RMSF plot with secondary structure annotation: `rmsf_plot.py`
2. RMSD calculation: `mda_rmsd.py`
3. Hydrogen bond calculation: `get_ave_hbond.py` 

## Ion channel analysis
1. Ion, water PMF calculation for ion channels:
    ```python
    ion_water_pmf.py
    ion_water_pmf_plt.py
    ```
2. CHAP (channel annotation) pore plot: `python chap_pore_plt.py`

## Secondary structure analysis
dir: `SecondaryStructureAnalysis`

1. DSSP in the `mdtraj` package: `mdtraj_dssp.py`
2. STRIDE in VMD: `vmd_stride.tcl`. Check: [Using VMD STRIDE](./SecondaryStructureAnalysis/README.md)
