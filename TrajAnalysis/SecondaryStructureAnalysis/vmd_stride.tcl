#!/usr/bin/tclsh
# VMD Secondary Structure Analysis Script
# Usage: vmd -dispdev text -e vmd_stride.tcl

# Configuration parameters - modify these as needed
set pdb_file "step4.1_equilibration.pdb"
set xtc_file "r1_500ns_303K_cyto_dt1ns.xtc"
set selection_text "chain B and resid 205 to 210"
set output_file "secondary_structure_analysis.txt"

# Load the structure and trajectory without graphics display
mol new $pdb_file type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set molid [molinfo top get id]

# Add trajectory
mol addfile $xtc_file type xtc first 0 last -1 step 10 filebonds 1 autobonds 1 waitfor all molid $molid

# Get number of frames
set num_frames [molinfo $molid get numframes]
puts "Loaded $num_frames frames for molecule ID: $molid"

# Test the selection first
set test_sel [atomselect $molid $selection_text]
set test_num [$test_sel num]
puts "Selection '$selection_text' contains $test_num atoms"
$test_sel delete

# Test CA selection
set test_sel_ca [atomselect $molid "($selection_text) and name CA"]
set test_num_ca [$test_sel_ca num]
puts "CA selection '($selection_text) and name CA' contains $test_num_ca atoms"

if {$test_num_ca == 0} {
    puts "ERROR: No CA atoms found! Check your selection."
    exit 1
}

$test_sel_ca delete

# Create selection
set sel [atomselect $molid $selection_text]
set residue_list [$sel get {chain resid}]

# Remove duplicates and sort residues
set unique_residues {}
foreach residue $residue_list {
    if {[lsearch $unique_residues $residue] == -1} {
        lappend unique_residues $residue
    }
}
set unique_residues [lsort -unique $unique_residues]

puts "Analyzing [llength $unique_residues] unique residues: $unique_residues"

# Initialize data structure to store results
array unset ss_data
foreach residue_pair $unique_residues {
    set chain [lindex $residue_pair 0]
    set resid [lindex $residue_pair 1]
    set key "${chain}:${resid}"
    set ss_data($key) {}
}

# Function to get secondary structure using VMD's built-in STRIDE
proc get_secondary_structure {molid frame selection_text} {
    # Set the frame
    animate goto $frame
    
    # Update secondary structure calculation using VMD's internal STRIDE
    mol ssrecalc $molid
    
    # Get the selection for this frame (select only CA atoms to avoid duplicates)
    set sel [atomselect $molid "($selection_text) and name CA" frame $frame]
    
    # Debug: Check if selection is valid
    set num_atoms [$sel num]
    if {$num_atoms == 0} {
        puts "Warning: No atoms found for selection '($selection_text) and name CA' in frame $frame"
        $sel delete
        return {}
    }
    
    # Get secondary structure data
    set ss_list [$sel get structure]
    set chain_list [$sel get chain]
    set resid_list [$sel get resid]
    
    # Debug: Print first few values for frame 0
    if {$frame == 0} {
        puts "Debug - Frame 0 data:"
        puts "  Number of CA atoms: $num_atoms"
        puts "  Chains: [lrange $chain_list 0 2]..."
        puts "  Resids: [lrange $resid_list 0 2]..."
        puts "  SS: [lrange $ss_list 0 2]..."
    }
    
    # Create a list to return results
    set result_list {}
    for {set i 0} {$i < [llength $chain_list]} {incr i} {
        set chain [lindex $chain_list $i]
        set resid [lindex $resid_list $i]
        set ss [lindex $ss_list $i]
        
        # Handle empty or undefined secondary structure
        if {$ss == "" || $ss == "0" || $ss == "{}"} {
            set ss "C"  ; # Default to coil if undefined
        }
        
        set key "${chain}:${resid}"
        lappend result_list $key $ss
    }
    
    $sel delete
    return $result_list
}

# Alternative function using Timeline plugin (if available)
proc get_secondary_structure_timeline {molid frame selection_text} {
    # This requires the Timeline plugin to be loaded
    # package require timeline
    
    animate goto $frame
    
    # Timeline plugin method would go here
    # For now, we'll use the standard VMD method
    return [get_secondary_structure $molid $frame $selection_text]
}

# Process each frame
puts "Processing frames..."
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {$frame % 100 == 0} {
        puts "Processing frame $frame/$num_frames"
    }
    
    # Get secondary structure for this frame
    set frame_results [get_secondary_structure $molid $frame $selection_text]
    
    # Debug: Print some results for first frame
    if {$frame == 0} {
        puts "Debug - Frame 0 results: [lrange $frame_results 0 5]..."
    }
    
    # Store results (frame_results is a flat list: key1 value1 key2 value2 ...)
    foreach {key ss_value} $frame_results {
        if {[info exists ss_data($key)]} {
            lappend ss_data($key) $ss_value
        } else {
            puts "Warning: Key $key not found in ss_data array"
        }
    }
    
    # Debug: Check data for first frame
    if {$frame == 0} {
        foreach key [lsort [array names ss_data]] {
            if {[llength $ss_data($key)] > 0} {
                puts "Debug - $key has data: [lindex $ss_data($key) 0]"
            } else {
                puts "Debug - $key has no data"
            }
        }
    }
}

# Write results to file
puts "Writing results to $output_file"
set output [open $output_file w]

# Write header - now residue names are columns
puts -nonewline $output "Frame"
foreach residue_pair $unique_residues {
    set chain [lindex $residue_pair 0]
    set resid [lindex $residue_pair 1]
    set residue_key "${chain}:${resid}"
    puts -nonewline $output "\t$residue_key"
}
puts $output ""

# Write data - one frame per row
for {set frame 0} {$frame < $num_frames} {incr frame} {
    puts -nonewline $output $frame
    
    foreach residue_pair $unique_residues {
        set chain [lindex $residue_pair 0]
        set resid [lindex $residue_pair 1]
        set residue_key "${chain}:${resid}"
        
        if {[info exists ss_data($residue_key)] && [llength $ss_data($residue_key)] > $frame} {
            set ss_value [lindex $ss_data($residue_key) $frame]
            puts -nonewline $output "\t$ss_value"
        } else {
            puts -nonewline $output "\t-"
        }
    }
    puts $output ""
}
exit
