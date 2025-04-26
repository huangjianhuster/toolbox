# Global variable to store the ID of the sphere
set sphere_id ""

# Set up a callback procedure to update the sphere position during animation
proc update_sphere_position {} {
    global sphere_id

    # global vmd_frame
    set selection [atomselect top "chain J and resname POPI and name C12 O2 C13 O3 C14 O4 P4 OP42 OP43 OP44 C15 O5 P5 OP52 OP53 OP54 C16 O6 C11 P O13 O14 O12 O11"]
    set com [measure center $selection weight mass]
    
    # Delete the old sphere
    if {$sphere_id ne ""} {
        draw delete $sphere_id
    }

    draw color green
    set sphere_id [draw sphere $com radius 1.0 resolution 50]

    # draw reference     
    # draw color blue
    # draw sphere {50 50 72} radius 2 resolution 50
    # draw sphere {50 90 72} radius 2 resolution 50
    # draw sphere {90 90 72} radius 2 resolution 50
    # draw sphere {86.2 92.7 67.9} radius 1.0 resolution 50
    # puts "drawing..."
}

# this function is required: https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node140.html
proc sphere_trace {args} {
    update_sphere_position
}

trace variable vmd_frame w sphere_trace
