proc default_rep {} {
    #### delete initial line representation for everything
    mol delrep 0 top
    
    #### make VDW representation for ions
    mol selection {ions}
    mol color Name
    mol representation VDW 
    mol material AOChalky
    mol addrep top

    mol selection {ions and ((x-70)**2+(y-70)**2) < 200}
    mol color Name
    mol representation VDW
    mol material AOChalky
    mol addrep top

    #### make Cartoon representation for protein
    mol selection {protein}
    mol color Chain
    mol representation NewCartoon
    mol material AOChalky
    mol addrep top

    ### add lipid
    mol selection {lipids and name P}
    mol color ColorID 2
    mol representation VDW 
    mol material AOChalky
    mol addrep top

    mol selection {lipids and not name P and not hydrogen}
    mol color ColorID 2
    mol representation Lines
    mol material AOChalky
    mol addrep top
    
    #### make line representation for water
    mol selection {water}
    mol color Name
    mol representation Lines
    mol material AOChalky
    mol addrep top
    
    set pbc_dims [pbc get -all]
    set pbc_size [lindex $pbc_dims 0]
    puts "PBC Box Size: $pbc_size"
    mol selection {same residue as (water and ((x-70)**2+(y-70)**2) < 200)}
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
