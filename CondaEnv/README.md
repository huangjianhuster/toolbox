conda/mamba environment configuration for computational biophysicists

1. install mamba
```bash
    # ref: https://github.com/conda-forge/miniforge
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh
    # following instructions
    # allow writing mamba into bashrc

    # of note, mamba from miniforge will automatic have the channel "conda-forge" prioritized
```

2. the "biosim" environment building process
```bash
    mamba create -n biosim python=3.10
    mamba activate biosim

    # unset
    unset $LD_LIBRARY_PATH

    # install basics
    mamba install numpy pandas matplotlib jupyterlab biopython mdanalysis mdtraj openmm biopython

    # install dssp in biopython
    mamba install -c salilab dssp

    # install pyrosetta
    # not for now...

    # psfgen
    # mamba install -c conda-forge psfgen

    # gnuplot
    mamba install conda-forge::gnuplot
    echo 'set term x11' > ~/.bashrc

    # of note, do not install vmd using mamba...it will interfere ipython

    # other setup
    # font in matplotlib
    sudo apt install msttcorefonts -qq
    rm ~/.cache/matplotlib -rf
```
 
