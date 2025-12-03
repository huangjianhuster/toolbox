# Global variable to store the ID of the sphere
set sphere_id [list]
set selection "chain A"

# Set up a callback procedure to update the sphere position during animation
proc update_sphere_position {} {
    global sphere_id
    global selection

    # global vmd_frame
    set selection [atomselect top $selection]
    set com [measure center $selection weight mass]
    
    # Delete the old sphere
    foreach id $sphere_id {
        draw delete $id
    }

    set sphere_id [list]

    draw color green
    set new_id [draw sphere $com radius 1.0 resolution 50]
    lappend sphere_id $new_id

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

# update_sphere_position
