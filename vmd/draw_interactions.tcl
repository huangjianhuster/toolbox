# VMD Interaction Drawing Script
# Draws interactions between residues with optional strength-based coloring

proc read_data {filename} {
    set data [list]
    set f [open $filename r]
    set line_num 0
    set has_strength -1  ;# -1 = unknown, 0 = no strength, 1 = has strength
    
    while {[gets $f line] != -1} {
        incr line_num
        
        # Skip empty lines
        set trimmed_line [string trim $line]
        if {$trimmed_line == ""} {
            continue
        }
        
        # Skip comment lines (lines starting with # or !)
        if {[string match "#*" $trimmed_line] || [string match "!*" $trimmed_line]} {
            continue
        }
        
        # Split by whitespace (handles both spaces and tabs, including multiple consecutive)
        set fields [regexp -all -inline {\S+} $line]
        
        # Determine number of columns
        set num_fields [llength $fields]
        
        # Check if we have at least 6 columns (chain1, resid1, resname1, chain2, resid2, resname2)
        if {$num_fields < 6} {
            puts "WARNING: Line $line_num has fewer than 6 columns, skipping: $line"
            puts "         Expected format: chainid_1 residue_index_1 residue_name_1 chainid_2 residue_index_2 residue_name_2 \[strength\]"
            continue
        }
        
        # Detect if we have strength column (first valid data line determines format)
        if {$has_strength == -1} {
            if {$num_fields >= 7} {
                # Check if 7th column is a valid number
                if {[string is double -strict [lindex $fields 6]]} {
                    set has_strength 1
                    puts "Detected format: chainid_1 resid_1 resname_1 chainid_2 resid_2 resname_2 interaction_strength"
                } else {
                    set has_strength 0
                    puts "Detected format: chainid_1 resid_1 resname_1 chainid_2 resid_2 resname_2 (no strength values)"
                }
            } else {
                set has_strength 0
                puts "Detected format: chainid_1 resid_1 resname_1 chainid_2 resid_2 resname_2 (no strength values)"
            }
        }
        
        # Parse residue 1 information
        set chain1 [lindex $fields 0]
        set resid1 [lindex $fields 1]
        set resname1 [lindex $fields 2]
        
        # Parse residue 2 information
        set chain2 [lindex $fields 3]
        set resid2 [lindex $fields 4]
        set resname2 [lindex $fields 5]
        
        # Validate residue IDs are numeric or alphanumeric
        if {![string is integer -strict $resid1] && ![regexp {^[0-9]+[A-Z]?$} $resid1]} {
            puts "WARNING: Line $line_num has invalid residue index for residue 1: $resid1"
            continue
        }
        
        if {![string is integer -strict $resid2] && ![regexp {^[0-9]+[A-Z]?$} $resid2]} {
            puts "WARNING: Line $line_num has invalid residue index for residue 2: $resid2"
            continue
        }
        
        # Get the strength value if present
        if {$has_strength == 1} {
            set strength [lindex $fields 6]
            
            # Validate that strength is a number
            if {![string is double -strict $strength]} {
                puts "WARNING: Line $line_num has non-numeric strength value: $strength, skipping"
                continue
            }
        } else {
            set strength -1  ;# Use -1 to indicate no strength value
        }
        
        # Store the parsed data
        lappend data [list $chain1 $resid1 $resname1 $chain2 $resid2 $resname2 $strength]
    }
    close $f
    
    puts "Successfully read [llength $data] interactions from $filename"
    if {$has_strength == 1} {
        puts "Interaction strengths will be used for coloring and scaling"
    } else {
        puts "No interaction strengths detected - uniform visualization will be used"
    }
    
    return $data
}

proc draw_interactions {data {default_color "blue"} {default_radius 0.2} {weak_cutoff 0.3} {strong_cutoff 0.7} {scale 1.0} {draw_node 0}} {
    
    # Color scheme for strength-based coloring
    # weak (< weak_cutoff): green
    # medium (weak_cutoff <= strength < strong_cutoff): orange
    # strong (>= strong_cutoff): red
    
    set interaction_count 0
    set failed_selections 0
    set has_strength_values 0
    set strength_min 1e10
    set strength_max -1e10
    
    # First pass: check if we have strength values and find min/max
    foreach interaction $data {
        set strength [lindex $interaction 6]
        if {$strength >= 0} {
            set has_strength_values 1
            if {$strength < $strength_min} {set strength_min $strength}
            if {$strength > $strength_max} {set strength_max $strength}
        }
    }
    
    if {$has_strength_values} {
        puts "\n========================================="
        puts "Drawing interactions with strength-based coloring"
        puts "Strength range: [format "%.3f" $strength_min] to [format "%.3f" $strength_max]"
        puts "Color scheme:"
        puts "  Weak   (< [format "%.2f" $weak_cutoff]): green"
        puts "  Medium ([format "%.2f" $weak_cutoff] - [format "%.2f" $strong_cutoff]): orange"
        puts "  Strong (>= [format "%.2f" $strong_cutoff]): red"
        puts "Radius scale factor: $scale"
        puts "========================================="
    } else {
        puts "\n========================================="
        puts "Drawing interactions with uniform style"
        puts "Color: $default_color"
        puts "Radius: [format "%.3f" $default_radius]"
        puts "========================================="
    }
    
    foreach interaction $data {
        set chain1 [lindex $interaction 0]
        set resid1 [lindex $interaction 1]
        set resname1 [lindex $interaction 2]
        set chain2 [lindex $interaction 3]
        set resid2 [lindex $interaction 4]
        set resname2 [lindex $interaction 5]
        set strength [lindex $interaction 6]
        
        # Create selection for first residue
        set sel1_string "chain $chain1 and resid $resid1 and name CA"
        set sel1 [atomselect top $sel1_string]
        
        if {[$sel1 num] == 0} {
            incr failed_selections
            puts "WARNING: Selection failed for residue 1 - No CA atom found"
            puts "         Chain: $chain1, ResID: $resid1, ResName: $resname1"
            puts "         Selection string: '$sel1_string'"
            puts "         Please verify:"
            puts "           - Chain ID '$chain1' exists in the structure"
            puts "           - Residue $resid1 exists in chain $chain1"
            puts "           - The residue has a CA atom"
            $sel1 delete
            continue
        }
        
        set coords1 [lindex [$sel1 get {x y z}] 0]
        set sel1x [lindex $coords1 0]
        set sel1y [lindex $coords1 1]
        set sel1z [lindex $coords1 2]
        
        # Verify the residue name matches (optional check)
        set actual_resname1 [lindex [$sel1 get resname] 0]
        if {$resname1 != "" && $actual_resname1 != $resname1} {
            puts "NOTE: Residue name mismatch for chain $chain1 resid $resid1"
            puts "      Expected: $resname1, Found: $actual_resname1"
        }
        
        # Create selection for second residue
        set sel2_string "chain $chain2 and resid $resid2 and name CA"
        set sel2 [atomselect top $sel2_string]
        
        if {[$sel2 num] == 0} {
            incr failed_selections
            puts "WARNING: Selection failed for residue 2 - No CA atom found"
            puts "         Chain: $chain2, ResID: $resid2, ResName: $resname2"
            puts "         Selection string: '$sel2_string'"
            puts "         Please verify:"
            puts "           - Chain ID '$chain2' exists in the structure"
            puts "           - Residue $resid2 exists in chain $chain2"
            puts "           - The residue has a CA atom"
            $sel1 delete
            $sel2 delete
            continue
        }
        
        set coords2 [lindex [$sel2 get {x y z}] 0]
        set sel2x [lindex $coords2 0]
        set sel2y [lindex $coords2 1]
        set sel2z [lindex $coords2 2]
        
        # Verify the residue name matches (optional check)
        set actual_resname2 [lindex [$sel2 get resname] 0]
        if {$resname2 != "" && $actual_resname2 != $resname2} {
            puts "NOTE: Residue name mismatch for chain $chain2 resid $resid2"
            puts "      Expected: $resname2, Found: $actual_resname2"
        }
        
        # Determine color and radius based on strength
        if {$has_strength_values && $strength >= 0} {
            # Strength-based coloring
            if {$strength < $weak_cutoff} {
                set color "green"
            } elseif {$strength < $strong_cutoff} {
                set color "orange"
            } else {
                set color "red"
            }
            # Use strength value as radius, scaled by the scale factor
            set radius [expr {$strength * $scale}]
        } else {
            # Uniform coloring when no strength values
            set color $default_color
            set radius $default_radius
        }
        
        # Draw the cylinder
        draw color $color
        draw cylinder "$sel1x $sel1y $sel1z" "$sel2x $sel2y $sel2z" radius $radius
        
        incr interaction_count

        # Clean up selections
        $sel1 delete
        $sel2 delete

        # Draw nodes if requested
        if {$draw_node == 1} {
            mol representation VDW 1.200000 0.000000
            mol color ColorID 1
            mol selection "chain $chain1 and resid $resid1 and name CA"
            mol material Opaque
            mol addrep top
            mol selection "chain $chain2 and resid $resid2 and name CA"
            mol material Opaque
            mol addrep top    
        }
    }
    
    puts "\n========================================="
    puts "Drawing Summary"
    puts "========================================="
    puts "Successfully drew: $interaction_count interactions"
    if {$failed_selections > 0} {
        puts "Failed selections: $failed_selections interactions"
        puts "\nSuggestions for failed selections:"
        puts "  1. Check that chain IDs in your data file match the loaded structure"
        puts "  2. Verify residue numbering matches between data and structure"
        puts "  3. Use VMD command 'set sel \[atomselect top \"chain X\"\]; \$sel get resid'"
        puts "     to see available residues in chain X"
        puts "  4. Check if the structure has CA atoms (some structures may only have backbone)"
    }
    puts "========================================="
}

# Helper function to list available chains and residues
proc list_structure_info {} {
    puts "\n========================================="
    puts "Structure Information"
    puts "========================================="
    
    set all [atomselect top "all"]
    set chains [lsort -unique [$all get chain]]
    
    puts "Available chains: $chains"
    puts ""
    
    foreach chain $chains {
        set chain_sel [atomselect top "chain $chain"]
        set resids [lsort -unique [$chain_sel get resid]]
        set num_residues [llength $resids]
        
        puts "Chain $chain:"
        puts "  Number of residues: $num_residues"
        if {$num_residues <= 20} {
            puts "  Residue IDs: $resids"
        } else {
            set first_10 [lrange $resids 0 9]
            set last_10 [lrange $resids end-9 end]
            puts "  First 10 residue IDs: $first_10"
            puts "  Last 10 residue IDs: $last_10"
        }
        
        # Check for CA atoms
        set ca_sel [atomselect top "chain $chain and name CA"]
        set num_ca [$ca_sel num]
        puts "  CA atoms: $num_ca"
        
        $chain_sel delete
        $ca_sel delete
        puts ""
    }
    
    $all delete
    puts "========================================="
}

# Print usage information
proc print_usage {} {
    puts "\n========================================="
    puts "VMD Interaction Drawing Script - Usage"
    puts "========================================="
    puts ""
    puts "Functions:"
    puts ""
    puts "1. read_data <filename>"
    puts "   Reads interaction data from a space/tab-separated file"
    puts ""
    puts "   Expected format (6 or 7 columns):"
    puts "     chainid_1 resid_1 resname_1 chainid_2 resid_2 resname_2 \[strength\]"
    puts ""
    puts "   - Lines starting with # or ! are treated as comments"
    puts "   - Empty lines are skipped"
    puts "   - Interaction strength (7th column) is optional"
    puts "   - Chain IDs and residue indices together provide specific selection"
    puts ""
    puts "2. draw_interactions <data> \[options\]"
    puts "   Draws cylinders between residue pairs"
    puts ""
    puts "   Parameters:"
    puts "     data            : output from read_data (required)"
    puts "     default_color   : color when no strength values (default: 'blue')"
    puts "     default_radius  : radius when no strength values (default: 0.2)"
    puts "     weak_cutoff     : cutoff for weak interactions (default: 0.3)"
    puts "     strong_cutoff   : cutoff for strong interactions (default: 0.7)"
    puts "     scale           : scaling factor for strength-based radius (default: 1.0)"
    puts "     draw_node       : draw VDW spheres at nodes (default: 0)"
    puts ""
    puts "   Color scheme (when strength values present):"
    puts "     Weak   (< weak_cutoff)                  : green"
    puts "     Medium (weak_cutoff - strong_cutoff)    : orange"
    puts "     Strong (>= strong_cutoff)               : red"
    puts ""
    puts "3. list_structure_info"
    puts "   Helper function to display available chains and residues"
    puts "   Useful for troubleshooting selection failures"
    puts ""
    puts "========================================="
    puts "Examples:"
    puts "========================================="
    puts ""
    puts "Example 1: Check structure information first"
    puts "  list_structure_info"
    puts ""
    puts "Example 2: Basic usage with strength values"
    puts "  set interactions \[read_data \"interactions.dat\"\]"
    puts "  draw_interactions \$interactions"
    puts ""
    puts "Example 3: Custom cutoffs for strength-based coloring"
    puts "  set interactions \[read_data \"interactions.dat\"\]"
    puts "  draw_interactions \$interactions blue 0.2 0.5 0.9 1.0 0"
    puts "  # weak < 0.5 (green), medium 0.5-0.9 (orange), strong >= 0.9 (red)"
    puts ""
    puts "Example 4: Without strength values (uniform)"
    puts "  set interactions \[read_data \"interactions_no_strength.dat\"\]"
    puts "  draw_interactions \$interactions red 0.3"
    puts "  # All interactions drawn in red with radius 0.3"
    puts ""
    puts "Example 5: With node visualization"
    puts "  set interactions \[read_data \"interactions.dat\"\]"
    puts "  draw_interactions \$interactions blue 0.2 0.3 0.7 1.5 1"
    puts "  # Scale radius by 1.5x and draw VDW spheres at nodes"
    puts ""
    puts "========================================="
    puts "Example data file formats:"
    puts "========================================="
    puts ""
    puts "Format 1: With interaction strengths"
    puts "  # Chain ResID ResName Chain ResID ResName Strength"
    puts "  A  10  ALA  B  25  GLY  0.25"
    puts "  A  15  LEU  B  30  VAL  0.55"
    puts "  C  45  SER  D  67  THR  0.85"
    puts ""
    puts "Format 2: Without interaction strengths"
    puts "  # Chain ResID ResName Chain ResID ResName"
    puts "  A  10  ALA  B  25  GLY"
    puts "  A  15  LEU  B  30  VAL"
    puts "  C  45  SER  D  67  THR"
    puts ""
    puts "Format 3: Mixed whitespace (tabs and spaces)"
    puts "  A	10	ALA	B	25	GLY	0.8"
    puts "  A  15  LEU  B  30  VAL  1.2"
    puts "  C    45    SER    D    67    THR    0.5"
    puts ""
    puts "Format 4: With insertion codes (e.g., 100A)"
    puts "  A  100A  ALA  B  25  GLY  0.75"
    puts "========================================="
}

puts "\n========================================="
puts "VMD Interaction Drawing Script Loaded"
puts "========================================="
puts "Useful commands:"
puts "  print_usage          - Show detailed usage information"
puts "  list_structure_info  - Show chains and residues in loaded structure"
puts "========================================="
