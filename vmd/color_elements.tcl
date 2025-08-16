proc trpv4_color {chainid1 chainid2} {
    #### delete initial line representation for everything
    # mol delrep 0 top

    mol selection "protein and segname $chainid1 $chainid2 and resid 146 to 366"
    mol color ColorID 6
    # mol representation NewCartoon [thickness] [resolution] [Aspect Ratio] [ID]
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 367 to 465 750 to 900"
    mol color ColorID 3
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 466 to 596"
    mol color ColorID 22
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 597 to 612"
    mol color ColorID 19
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 613 to 638"
    mol color ColorID 5
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 662 to 706"
    mol color ColorID 13
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 707 to 724"
    mol color ColorID 27
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 725 to 735"
    mol color ColorID 0
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top
    mol selection "protein and segname $chainid1 $chainid2 and resid 736 to 749"
    mol color ColorID 30
    mol representation NewCartoon 0.300000 10.000000 3.080000 0
    mol material Fading
    mol addrep top

    #### turn off axes
    axes location off

    ### set background white
    # color Display Background 8

}

proc clean_rep {} {
    set numReps [molinfo top get numreps]
    for {set i [expr {$numReps - 1}]} {$i >= 0} {incr i -1} {
        mol delrep $i top
    }
}


proc change_cartoon_material {new_material} {
    set molID [molinfo top]  ;# Get the top molecule ID
    puts $molID
    if {$molID < 0} {
        puts "No molecule is loaded."
        return
    }

    set numReps [molinfo $molID get numreps]
    # puts $numReps 
    for {set i 0} {$i < $numReps} {incr i} {
        
        set drawMethod [lindex [mol modstyle $molID $i] 0]
        puts $drawMethod 
        if {[string match "NewCartoon" $drawMethod]} {
            mol modmaterial $molID $i $new_material
        }
    }
}

