proc default_rep {} {
    #### delete initial line representation for everything
    mol delrep 0 top
    
    #### make VDW representation for ions
    mol selection {ions}
    mol color Name
    mol representation VDW 
    mol material AOChalky
    mol addrep top
    
    #### make Cartoon representation for protein
    mol selection {protein}
    mol color Structure
    mol representation NewCartoon
    mol material AOChalky
    mol addrep top
    
    #### make line representation for water
    mol selection {water}
    mol color Name
    mol representation Lines
    mol material AOChalky
    mol addrep top
    
    #### turn off axes
    axes location off
    
    ### set background white
    color Display Background 8

    ### pbc box
    pbc box -width 2.0
}

default_rep
