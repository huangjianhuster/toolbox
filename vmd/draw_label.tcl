# First, define a global list to store graphics IDs (to delete them later when tracing frames)
set text_labels {}

# Three letter to One letter mapping
set aa1letter {
    ALA A    ARG R    ASN N    ASP D    CYS C    GLN Q    GLU E    GLY G
    HIS H    HIE H    HSE H    HID H    ILE I    LEU L    LYS K    MET M
    PHE F    PRO P    SER S    THR T    TRP W    TYR Y    VAL V    HSD H
}

# Selection definitions: better to keep the variable name as "reslabel"; in the selection, need to include "name CA"
set reslabel [atomselect top "resid 450 to 750 and name CA"]


# Define a procedure to update the labels
proc update_labels {} {
    global text_labels
    global aa1letter
    global reslabel

    # Delete old labels
    foreach id $text_labels {
        graphics top delete $id
    }
    set text_labels {}

    # Select contacting residues (adjust your ligand selection if needed)
    $reslabel update

    # Get list of contacting residues
    set reslist [$reslabel get {resname resid x y z}]
    # puts $reslist

    # Define groups
    set negative {ASP GLU}
    set positive {LYS ARG}
    set polar {SER THR ASN GLN TYR CYS HIS HIE HSE}
    set hydrophobic {ALA VAL LEU ILE MET PHE TRP PRO GLY}
    
    # Loop through all residues and labeling them
    foreach res $reslist {

        lassign $res resname resid x y z

        # label shift
        set x_shift [expr {$x + 1}]
        set y_shift [expr {$y + 1}]
        set z_shift [expr {$z + 1}]

        # Determine color based on residue type
        # if {[lsearch $negative $resname] >= 0} {
        #     set labelcolor red
        # } elseif {[lsearch $positive $resname] >= 0} {
        #     set labelcolor blue
        # } elseif {[lsearch $polar $resname] >= 0} {
        #     set labelcolor cyan
        # } elseif {[lsearch $hydrophobic $resname] >= 0} {
        #     set labelcolor white
        # } else {
        #     set labelcolor yellow  ;# fallback color for unknown residues
        # }

        # override colors for labels if you want
        # Set label color blue (ColorID 0 = blue by default)
        # graphics top color 0
        set labelcolor black

        draw color $labelcolor

        # map to one letter
        set clean_resname [string trim $resname]

        if {[dict exists $aa1letter $clean_resname]} {
            set oneletter [dict get $aa1letter $clean_resname]
        } else {
            set oneletter "X"  ;# fallback for unknown residue
        }

        set label "${oneletter}${resid}"

        # Three letter version
        # set label "${resname}${resid}"

        # Draw the text
        set id [graphics top text "$x_shift $y_shift $z_shift" $label size 0.7 thickness 2.5]
        lappend text_labels $id
    }
}

proc label_trace {args} {
    update_labels
}

proc remove_labels {} {
    graphics top delete all
}

# Attach the update to frame changes
trace variable vmd_frame w label_trace
