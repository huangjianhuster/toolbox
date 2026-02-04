#!/usr/bin/tclsh
# VMD Alignment Script with Ligand Extraction
# Usage: source alignment_script.tcl
#        align_and_save_ligand [options]

proc align_and_save_ligand {args} {
    # Default parameters
    set mol_move 0
    set mol_target 1
    set align_sel "name CA"
    set align_sel_move ""
    set align_sel_target ""
    set ligand_sel "resname LIG"
    set output_file "aligned_ligand.pdb"
    set reverse 0
    
    # Parse arguments
    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]
        switch -- $arg {
            "-align" {
                incr i
                set align_sel [lindex $args $i]
            }
            "-align_move" {
                incr i
                set align_sel_move [lindex $args $i]
            }
            "-align_target" {
                incr i
                set align_sel_target [lindex $args $i]
            }
            "-ligand" {
                incr i
                set ligand_sel [lindex $args $i]
            }
            "-output" {
                incr i
                set output_file [lindex $args $i]
            }
            "-reverse" {
                set reverse 1
            }
            "-mol_move" {
                incr i
                set mol_move [lindex $args $i]
            }
            "-mol_target" {
                incr i
                set mol_target [lindex $args $i]
            }
            "-help" {
                print_usage
                return
            }
            default {
                puts "Unknown option: $arg"
                print_usage
                return
            }
        }
    }
    
    # Determine final selection strings
    # If separate selections are provided, use them; otherwise use the common selection
    if {$align_sel_move == ""} {
        set align_sel_move $align_sel
    }
    if {$align_sel_target == ""} {
        set align_sel_target $align_sel
    }
    
    # Handle reverse option
    if {$reverse} {
        set temp $mol_move
        set mol_move $mol_target
        set mol_target $temp
        
        set temp_sel $align_sel_move
        set align_sel_move $align_sel_target
        set align_sel_target $temp_sel
        
        puts "Reverse mode: Moving molecule $mol_move to align with molecule $mol_target"
    }
    
    puts "========================================="
    puts "VMD Alignment Script"
    puts "========================================="
    puts "Moving molecule: $mol_move"
    puts "Target molecule: $mol_target"
    puts "Alignment selection (moving): '$align_sel_move'"
    puts "Alignment selection (target): '$align_sel_target'"
    puts "Ligand selection: '$ligand_sel'"
    puts "Output file: $output_file"
    puts "========================================="
    
    # Check if molecules exist
    if {[molinfo num] < 2} {
        puts "ERROR: Need at least 2 molecules loaded!"
        puts "Current number of molecules: [molinfo num]"
        return -1
    }
    
    if {$mol_move >= [molinfo num] || $mol_target >= [molinfo num]} {
        puts "ERROR: Invalid molecule ID!"
        puts "Valid molecule IDs: 0 to [expr [molinfo num] - 1]"
        return -1
    }
    
    # Create selections for alignment
    puts "\n--- Creating alignment selections ---"
    
    if {[catch {set sel_target [atomselect $mol_target $align_sel_target]} err]} {
        puts "ERROR: Failed to create target selection with '$align_sel_target'"
        puts "Error message: $err"
        return -1
    }
    
    if {[catch {set sel_move [atomselect $mol_move $align_sel_move]} err]} {
        puts "ERROR: Failed to create moving selection with '$align_sel_move'"
        puts "Error message: $err"
        $sel_target delete
        return -1
    }
    
    set num_target [$sel_target num]
    set num_move [$sel_move num]
    
    puts "Target selection (mol $mol_target, '$align_sel_target'): $num_target atoms"
    puts "Moving selection (mol $mol_move, '$align_sel_move'): $num_move atoms"
    
    # Check if selections are valid
    if {$num_target == 0} {
        puts "ERROR: Target selection is empty!"
        puts "Please verify your selection string: '$align_sel_target'"
        $sel_target delete
        $sel_move delete
        return -1
    }
    
    if {$num_move == 0} {
        puts "ERROR: Moving selection is empty!"
        puts "Please verify your selection string: '$align_sel_move'"
        $sel_target delete
        $sel_move delete
        return -1
    }
    
    # Check if selections have the same number of atoms
    if {$num_target != $num_move} {
        puts "ERROR: Selection atom counts don't match!"
        puts "  Target (mol $mol_target, '$align_sel_target'): $num_target atoms"
        puts "  Moving (mol $mol_move, '$align_sel_move'): $num_move atoms"
        puts "VMD requires equal atom counts for alignment."
        puts ""
        puts "Suggestions:"
        puts "  - Check that both selections use compatible atom filters"
        puts "  - Verify residue numbering/naming schemes match"
        puts "  - Try using 'name CA' for both if working with proteins"
        $sel_target delete
        $sel_move delete
        return -1
    }
    
    puts "✓ Selections are valid and compatible"
    
    # Display sample atoms for verification
    puts "\n--- Sample atoms from selections ---"
    if {$num_target > 0} {
        set sample_target_atoms [lrange [$sel_target get {segname resname resid name}] 0 4]
        puts "Target selection (first 5 atoms):"
        set atom_idx 0
        foreach atom_info $sample_target_atoms {
            lassign $atom_info segname resname resid name
            puts "  Atom [expr $atom_idx + 1]: $segname $resname $resid $name"
            incr atom_idx
        }
        if {$num_target > 5} {
            puts "  ... and [expr $num_target - 5] more atoms"
        }
    }
    
    if {$num_move > 0} {
        set sample_move_atoms [lrange [$sel_move get {segname resname resid name}] 0 4]
        puts "\nMoving selection (first 5 atoms):"
        set atom_idx 0
        foreach atom_info $sample_move_atoms {
            lassign $atom_info segname resname resid name
            puts "  Atom [expr $atom_idx + 1]: $segname $resname $resid $name"
            incr atom_idx
        }
        if {$num_move > 5} {
            puts "  ... and [expr $num_move - 5] more atoms"
        }
    }
    
    # Perform alignment
    puts "\n--- Performing alignment ---"
    set sel_move_whole [atomselect $mol_move "all"]
    
    if {[catch {set transformation_matrix [measure fit $sel_move $sel_target]} err]} {
        puts "ERROR: Alignment failed!"
        puts "Error message: $err"
        $sel_target delete
        $sel_move delete
        $sel_move_whole delete
        return -1
    }
    
    $sel_move_whole move $transformation_matrix
    
    set rmsd [measure rmsd $sel_move $sel_target]
    puts "✓ Alignment complete"
    puts "  RMSD after alignment: [format "%.3f" $rmsd] Å"
    
    # Extract and save ligand
    puts "\n--- Extracting ligand ---"
    
    if {[catch {set sel_ligand [atomselect $mol_move $ligand_sel]} err]} {
        puts "ERROR: Failed to create ligand selection with '$ligand_sel'"
        puts "Error message: $err"
        $sel_target delete
        $sel_move delete
        $sel_move_whole delete
        return -1
    }
    
    set num_ligand [$sel_ligand num]
    puts "Ligand selection (mol $mol_move, '$ligand_sel'): $num_ligand atoms"
    
    if {$num_ligand == 0} {
        puts "WARNING: Ligand selection is empty!"
        puts "Please verify your ligand selection string: '$ligand_sel'"
        puts "No ligand file will be saved."
    } else {
        puts "✓ Ligand selection is valid"
        
        # Display ligand information
        set ligand_resnames [lsort -unique [$sel_ligand get resname]]
        set ligand_resids [lsort -unique [$sel_ligand get resid]]
        puts "  Residue names: $ligand_resnames"
        puts "  Residue IDs: $ligand_resids"
        
        # Save ligand coordinates
        if {[catch {$sel_ligand writepdb $output_file} err]} {
            puts "ERROR: Failed to write ligand PDB file"
            puts "Error message: $err"
        } else {
            puts "✓ Ligand coordinates saved to: $output_file"
        }
    }
    
    # Cleanup
    $sel_target delete
    $sel_move delete
    $sel_move_whole delete
    $sel_ligand delete
    
    puts "\n========================================="
    puts "Alignment completed successfully!"
    puts "========================================="
    
    return 0
}

proc print_usage {} {
    puts "\n========================================="
    puts "VMD Alignment Script - Usage"
    puts "========================================="
    puts "align_and_save_ligand \[options\]"
    puts ""
    puts "Options:"
    puts "  -align <selection>         Common alignment selection for both molecules"
    puts "                             (default: 'name CA')"
    puts "  -align_move <selection>    Alignment selection for moving molecule"
    puts "                             (overrides -align for moving molecule)"
    puts "  -align_target <selection>  Alignment selection for target molecule"
    puts "                             (overrides -align for target molecule)"
    puts "  -ligand <selection>        Ligand selection (default: 'resname LIG')"
    puts "  -output <filename>         Output PDB file (default: 'aligned_ligand.pdb')"
    puts "  -reverse                   Reverse alignment (align B to A instead of A to B)"
    puts "  -mol_move <id>             Molecule ID to move (default: 0)"
    puts "  -mol_target <id>           Target molecule ID (default: 1)"
    puts "  -help                      Show this help message"
    puts ""
    puts "Note: If -align_move and -align_target are both specified, -align is ignored."
    puts "      This allows different selection strings for each molecule while"
    puts "      ensuring equal atom counts for proper alignment."
    puts ""
    puts "Examples:"
    puts ""
    puts "1. Basic usage (same selection for both molecules):"
    puts "   align_and_save_ligand"
    puts "   align_and_save_ligand -align \"protein and backbone\""
    puts ""
    puts "2. Different selections for each molecule:"
    puts "   align_and_save_ligand -align_move \"chain A and name CA\" \\"
    puts "                         -align_target \"chain B and name CA\""
    puts ""
    puts "3. Align specific residue ranges (different numbering):"
    puts "   align_and_save_ligand -align_move \"resid 10 to 100 and name CA\" \\"
    puts "                         -align_target \"resid 50 to 140 and name CA\""
    puts ""
    puts "4. Custom ligand and output:"
    puts "   align_and_save_ligand -ligand \"resname ATP\" \\"
    puts "                         -output \"ATP_aligned.pdb\""
    puts ""
    puts "5. Reverse alignment with different chain selections:"
    puts "   align_and_save_ligand -reverse \\"
    puts "                         -align_move \"chain X and name CA\" \\"
    puts "                         -align_target \"chain Y and name CA\""
    puts ""
    puts "6. Align backbone atoms excluding terminal residues:"
    puts "   align_and_save_ligand -align_move \"resid 5 to 95 and backbone\" \\"
    puts "                         -align_target \"resid 15 to 105 and backbone\""
    puts "========================================="
}

# Print usage information when script is loaded
puts "\n========================================="
puts "VMD Alignment Script Loaded"
puts "========================================="
puts "To run alignment with default settings:"
puts "  align_and_save_ligand"
puts ""
puts "For help and options:"
puts "  align_and_save_ligand -help"
puts "========================================="
