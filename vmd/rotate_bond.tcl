proc rotate_bond_C1star_N9 {molID angle_deg} {
    # rotate bond between C1* and N9
    set sel_rot [atomselect $molID "name N9 C8 H80 N7 C5 C6 N6 H61 H60 N1 C2 H2 N3 C4" ]
    set at_fix [atomselect $molID "name 'C1*'"]
    
    # get rotation axis and origin
    set sel_fixed [atomselect $molID "name 'C1*'"]
    set sel_anchor [atomselect $molID "name N9"]
    set p1 [lindex [$sel_fixed get {x y z}] 0]
    set p2 [lindex [$sel_anchor get {x y z}] 0]
    
    $sel_rot move [trans bond $p1 $p2 $angle_deg deg]

    #  clean up selections
    $sel_rot delete
    $sel_fixed delete
    $sel_anchor delete
}
