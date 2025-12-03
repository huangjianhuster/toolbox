set draw_id [list]

# define functions
proc draw_principal_axes {selection_string {scale 20} {arrow_radius 0.8}} {
    global draw_id

    # Create atom selection
    set sel [atomselect top $selection_string]
    
    # Check if selection is valid
    if {[$sel num] == 0} {
        puts "Error: No atoms selected with '$selection_string'"
        $sel delete
        return
    }

    # Delete the old drawings
    foreach id $draw_id {
        draw delete $id
    }
    set draw_id [list]

    # Get center of mass
    set com [measure center $sel weight mass]
    
    # Get inertia tensor and principal axes
    set inertia [measure inertia $sel]
    set axes [lindex $inertia 1]
    
    # Draw shortest axis (blue)
    set axis1 [lindex $axes 0]
    set vscaled [vecscale $scale $axis1]
    set end [vecadd $com $vscaled]
    draw color blue
    set arrow_ids [vmd_draw_arrow 0 $com $end $arrow_radius]
    set draw_id [concat $draw_id $arrow_ids]
    
    # Draw second longest axis (green)
    set axis2 [lindex $axes 1]
    set vscaled [vecscale $scale $axis2]
    set end [vecadd $com $vscaled]
    draw color green
    set arrow_ids [vmd_draw_arrow 0 $com $end $arrow_radius]
    set draw_id [concat $draw_id $arrow_ids]
    
    # Draw longest axis (red)
    set axis3 [lindex $axes 2]
    set vscaled [vecscale $scale $axis3]
    set end [vecadd $com $vscaled]
    draw color red
    set arrow_ids [vmd_draw_arrow 0 $com $end $arrow_radius]
    set draw_id [concat $draw_id $arrow_ids]
    
    # Clean up
    $sel delete
    
    puts "Principal axes drawn for selection: $selection_string"
    puts "Center of mass: $com"
}

# Keep your vmd_draw_arrow function
proc vmd_draw_arrow {mol start end {cradius 0.15}} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    set cylinder_id [graphics $mol cylinder $start $middle radius $cradius resolution 50]
    set cone_radius [expr {2 * $cradius}]
    set cone_id [graphics $mol cone $middle $end radius $cone_radius resolution 50]
    
    return [list $cylinder_id $cone_id]
}

set selection "chain A and backbone"

# draw arrow
draw_principal_axes $selection 20 0.8
