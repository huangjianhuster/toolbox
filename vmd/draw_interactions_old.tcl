
# Define the specific color to use
# set base_color "red"

proc read_data {filename} {
    set data [list]
    set f [open $filename r]
    while {[gets $f line] != -1} {
        set fields [split $line "\t"]
        set res1 [split [lindex $fields 0] ":"]
        set res2 [split [lindex $fields 1] ":"]
        lappend data [list [lindex $res1 0] [lindex $res1 1] [lindex $res1 2] [lindex $res2 0] [lindex $res2 1] [lindex $res2 2] [lindex $fields 2]]
    }
    close $f
    return $data
}

proc draw_interactions {data color {scale 1} {draw_node 0} } {
    
    # Define cutoffs for different shades of a specific color
    # set shade_cutoffs {0.3 0.6 0.9}
    # set base_color "red"
     
    foreach interaction $data {
        set chain1 [lindex $interaction 0]
        set resname1 [lindex $interaction 1]
        set resid1 [lindex $interaction 2]
        set chain2 [lindex $interaction 3]
        set resname2 [lindex $interaction 4]
        set resid2 [lindex $interaction 5]
        set radius [lindex $interaction 6]
        set strength [lindex $interaction 6]
        
        set sel1 [atomselect top "chain $chain1 and resid $resid1 and name CA"]
        set sel1x [$sel1 get x]
        set sel1y [$sel1 get y]
        set sel1z [$sel1 get z]
        set sel2 [atomselect top "chain $chain2 and resid $resid2 and name CA"]
        set sel2x [$sel2 get x]
        set sel2y [$sel2 get y]
        set sel2z [$sel2 get z]
        # Determine the shade of the color based on the interaction strength value
        # set shade "0.5"
        # foreach cutoff $shade_cutoffs s {0.6 0.8 1.0} {
        #     if {$strength >= $cutoff} {
        #         set shade $s
        #     } else {
        #         break
        #     }
        #     }

        # Construct the color string
        # set color "${base_color}${shade}"
      
        draw color $color
        set radius [expr {$radius * $scale}]
        # Create the cylinder using the residue numbers, radius, and color
        draw cylinder "$sel1x $sel1y $sel1z" "$sel2x $sel2y $sel2z" radius $radius

        # if draw_node == 1
        if {$draw_node == 1} {
            mol representation VDW 1.200000 0.000000
            mol color ColorID 1
            mol selection "chain $chain1 and resid $resid1 and name CA"
            # mol selection $sel2
            mol material Opaque
            mol addrep top
            mol selection "chain $chain2 and resid $resid2 and name CA"
            mol material Opaque
            mol addrep top    
        }

    }
}
