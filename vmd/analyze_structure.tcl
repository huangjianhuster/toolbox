# VMD Structure Analysis Script
# Provides detailed information about the loaded molecular structure

proc analyze_structure {{molid "top"}} {
    
    # Check if molecule exists
    if {$molid == "top"} {
        set molid [molinfo top]
    }
    
    if {$molid == -1} {
        puts "ERROR: No molecule loaded!"
        return
    }
    
    puts "\n========================================================================"
    puts "STRUCTURE ANALYSIS"
    puts "========================================================================"
    puts "Molecule ID: $molid"
    puts "Molecule name: [molinfo $molid get name]"
    
    set all [atomselect $molid "all"]
    set total_atoms [$all num]
    puts "Total atoms: $total_atoms"
    $all delete
    
    puts "\n------------------------------------------------------------------------"
    puts "CHAIN INFORMATION"
    puts "------------------------------------------------------------------------"
    
    # Get all chains
    set all_sel [atomselect $molid "all"]
    set all_chains [lsort -unique [$all_sel get chain]]
    $all_sel delete
    
    if {[llength $all_chains] == 0 || ([llength $all_chains] == 1 && [lindex $all_chains 0] == "")} {
        puts "No chain information available (single chain or no chain IDs)"
        set all_chains [list ""]
    }
    
    puts "Number of chains: [llength $all_chains]"
    if {[llength $all_chains] <= 50} {
        puts "Chain IDs: $all_chains"
    } else {
        puts "Chain IDs: Too many to display ([llength $all_chains] chains)"
    }
    puts ""
    
    # Analyze each chain
    foreach chain $all_chains {
        if {$chain == ""} {
            set chain_sel [atomselect $molid "all"]
            set chain_label "No chain ID"
        } else {
            set chain_sel [atomselect $molid "chain $chain"]
            set chain_label "Chain $chain"
        }
        
        set num_atoms [$chain_sel num]
        
        if {$num_atoms > 0} {
            puts "  $chain_label:"
            puts "    Total atoms: $num_atoms"
            
            # Get residue information
            set resids [lsort -unique [$chain_sel get resid]]
            set num_residues [llength $resids]
            puts "    Number of residues: $num_residues"
            
            if {$num_residues > 0} {
                # Find min and max resid
                set min_resid [lindex $resids 0]
                set max_resid [lindex $resids end]
                puts "    Residue index range: $min_resid to $max_resid"
                
                # Display residue list (smart display based on number)
                if {$num_residues <= 20} {
                    puts "    Residue indices: $resids"
                } else {
                    set first_10 [lrange $resids 0 9]
                    set last_10 [lrange $resids end-9 end]
                    puts "    First 10 residue indices: $first_10"
                    puts "    Last 10 residue indices: $last_10"
                }
            }
            
            # Check for protein
            set protein_sel [atomselect $molid "chain $chain and protein"]
            set num_protein_atoms [$protein_sel num]
            if {$num_protein_atoms > 0} {
                set protein_resids [lsort -unique [$protein_sel get resid]]
                puts "    Protein residues: [llength $protein_resids]"
            }
            $protein_sel delete
            
            # Check for nucleic acids
            set nucleic_sel [atomselect $molid "chain $chain and nucleic"]
            set num_nucleic_atoms [$nucleic_sel num]
            if {$num_nucleic_atoms > 0} {
                set nucleic_resids [lsort -unique [$nucleic_sel get resid]]
                puts "    Nucleic acid residues: [llength $nucleic_resids]"
            }
            $nucleic_sel delete
            
            puts ""
        }
        
        $chain_sel delete
    }
    
    puts "------------------------------------------------------------------------"
    puts "WATER MOLECULES"
    puts "------------------------------------------------------------------------"
    
    # Common water residue names
    set water_names [list "WAT" "HOH" "TIP3" "TIP4" "TIP5" "SPC" "SPCE" "H2O" "SOL"]
    set total_water_atoms 0
    set found_water_types [list]
    
    foreach water_name $water_names {
        set water_sel [atomselect $molid "resname $water_name"]
        set num_water_atoms [$water_sel num]
        if {$num_water_atoms > 0} {
            set num_water_molecules [expr {$num_water_atoms / 3}]  ;# Assuming 3 atoms per water
            puts "  $water_name: $num_water_molecules molecules ($num_water_atoms atoms)"
            lappend found_water_types $water_name
            set total_water_atoms [expr {$total_water_atoms + $num_water_atoms}]
        }
        $water_sel delete
    }
    
    if {$total_water_atoms == 0} {
        puts "  No water molecules detected"
    } else {
        puts "  ---"
        puts "  Total water atoms: $total_water_atoms (approximately [expr {$total_water_atoms / 3}] molecules)"
    }
    
    puts "\n------------------------------------------------------------------------"
    puts "IONS"
    puts "------------------------------------------------------------------------"
    
    # Common ion names
    set ion_types [list \
        [list "NA" "'Na+'" "SOD"] \
        [list "K" "'K+'" "POT"] \
        [list "CL" "'Cl-'" "CLA"] \
        [list "CA" "'Ca2+'" "CAL"] \
        [list "MG" "'Mg2+'"] \
        [list "ZN" "'Zn2+'"] \
        [list "FE" "'Fe2+'" "'Fe3+'"] \
        [list "CU" "'Cu+'" "'Cu2+'"] \
        [list "MN" "'Mn2+'"] \
    ]
    
    set total_ions 0
    set found_ions [dict create]
    
    foreach ion_group $ion_types {
        set ion_label [lindex $ion_group end]
        set ion_names [lrange $ion_group 0 end-1]
        set ion_count 0
        set found_resnames [list]
        
        foreach ion_name $ion_names {
            set ion_sel [atomselect $molid "resname $ion_name"]
            set num_ions [$ion_sel num]
            if {$num_ions > 0} {
                set ion_count [expr {$ion_count + $num_ions}]
                lappend found_resnames $ion_name
            }
            $ion_sel delete
        }
        
        if {$ion_count > 0} {
            puts "  $ion_label ([join $found_resnames {, }]): $ion_count ions"
            set total_ions [expr {$total_ions + $ion_count}]
        }
    }
    
    if {$total_ions == 0} {
        puts "  No ions detected"
    } else {
        puts "  ---"
        puts "  Total ions: $total_ions"
    }
    
    puts "\n------------------------------------------------------------------------"
    puts "LIGANDS AND SMALL MOLECULES"
    puts "------------------------------------------------------------------------"
    
    # Select everything that is not protein, nucleic, water, or ions
    set ligand_sel [atomselect $molid "not protein and not nucleic and not water and not ions"]
    set num_ligand_atoms [$ligand_sel num]
    
    if {$num_ligand_atoms > 0} {
        # Get unique residue names for ligands
        set ligand_resnames [lsort -unique [$ligand_sel get resname]]
        
        puts "  Detected non-standard residues (potential ligands/cofactors):"
        puts ""
        
        # Count each ligand type
        foreach resname $ligand_resnames {
            set this_ligand_sel [atomselect $molid "resname $resname"]
            set num_atoms [$this_ligand_sel num]
            set resids [lsort -unique [$this_ligand_sel get resid]]
            set num_molecules [llength $resids]
            
            # Get chains where this ligand is found
            set ligand_chains [lsort -unique [$this_ligand_sel get chain]]
            set chain_info [join $ligand_chains ", "]
            if {$chain_info == ""} {
                set chain_info "no chain ID"
            }
            
            puts "    $resname:"
            puts "      Molecules: $num_molecules"
            puts "      Total atoms: $num_atoms"
            puts "      Atoms per molecule: [expr {$num_atoms / $num_molecules}]"
            puts "      Found in chains: $chain_info"
            
            $this_ligand_sel delete
        }
        
        puts "  ---"
        puts "  Total ligand/small molecule atoms: $num_ligand_atoms"
        puts "  Number of different ligand types: [llength $ligand_resnames]"
    } else {
        puts "  No ligands or small molecules detected"
    }
    
    $ligand_sel delete
    
    puts "\n------------------------------------------------------------------------"
    puts "LIPIDS"
    puts "------------------------------------------------------------------------"
    
    # Common lipid residue names (CHARMM, AMBER, GROMOS naming conventions)
    set lipid_names [list \
        "POPC" "POPE" "POPS" "POPG" "DPPC" "DPPE" "DMPC" "DOPC" "DOPE" \
        "DLPC" "DLPE" "DSPC" "PALM" "STEA" "OLEO" \
        "CHOL" "CHL1" "CLOL" \
        "LPPC" "PAPC" "DAPC" "PGPC" "SAPC" \
        "PC" "PE" "PS" "PG" "PA" "PI" "SM" "CER" \
    ]
    
    set total_lipid_atoms 0
    set found_lipids [dict create]
    
    puts "  Checking for common lipid types..."
    puts ""
    
    foreach lipid_name $lipid_names {
        set lipid_sel [atomselect $molid "resname $lipid_name"]
        set num_lipid_atoms [$lipid_sel num]
        if {$num_lipid_atoms > 0} {
            set resids [lsort -unique [$lipid_sel get resid]]
            set num_lipid_molecules [llength $resids]
            puts "    $lipid_name: $num_lipid_molecules molecules ($num_lipid_atoms atoms)"
            set total_lipid_atoms [expr {$total_lipid_atoms + $num_lipid_atoms}]
            dict set found_lipids $lipid_name $num_lipid_molecules
        }
        $lipid_sel delete
    }
    
    if {$total_lipid_atoms == 0} {
        puts "  No lipids detected"
    } else {
        puts "  ---"
        puts "  Total lipid atoms: $total_lipid_atoms"
        puts "  Total lipid molecules: [expr [join [dict values $found_lipids] +]]"
    }
    
    puts "\n------------------------------------------------------------------------"
    puts "SECONDARY STRUCTURE (if assigned)"
    puts "------------------------------------------------------------------------"
    
    set protein_sel [atomselect $molid "protein"]
    if {[$protein_sel num] > 0} {
        # Check if secondary structure is assigned
        set ss_types [lsort -unique [$protein_sel get structure]]
        
        if {[llength $ss_types] > 1 || ([llength $ss_types] == 1 && [lindex $ss_types 0] != "C")} {
            # Count different secondary structures
            set helix_sel [atomselect $molid "protein and structure H"]
            set sheet_sel [atomselect $molid "protein and structure E"]
            set turn_sel [atomselect $molid "protein and structure T"]
            set coil_sel [atomselect $molid "protein and structure C"]
            
            puts "  Alpha helix (H): [$helix_sel num] atoms"
            puts "  Beta sheet (E): [$sheet_sel num] atoms"
            puts "  Turn (T): [$turn_sel num] atoms"
            puts "  Coil (C): [$coil_sel num] atoms"
            
            $helix_sel delete
            $sheet_sel delete
            $turn_sel delete
            $coil_sel delete
        } else {
            puts "  Secondary structure not assigned"
            puts "  To assign: mol ssrecalc top (in VMD console)"
        }
    } else {
        puts "  No protein in structure"
    }
    $protein_sel delete
    
    puts "\n------------------------------------------------------------------------"
    puts "SUMMARY"
    puts "------------------------------------------------------------------------"
    
    # Overall composition
    set protein_all [atomselect $molid "protein"]
    set nucleic_all [atomselect $molid "nucleic"]
    set water_all [atomselect $molid "water"]
    set ion_all [atomselect $molid "ions"]
    
    puts "  Total atoms: $total_atoms"
    puts "    Protein: [$protein_all num] atoms ([format "%.1f" [expr {100.0 * [$protein_all num] / $total_atoms}]]%)"
    puts "    Nucleic acid: [$nucleic_all num] atoms ([format "%.1f" [expr {100.0 * [$nucleic_all num] / $total_atoms}]]%)"
    puts "    Water: [$water_all num] atoms ([format "%.1f" [expr {100.0 * [$water_all num] / $total_atoms}]]%)"
    puts "    Ions: [$ion_all num] atoms ([format "%.1f" [expr {100.0 * [$ion_all num] / $total_atoms}]]%)"
    
    set other_atoms [expr {$total_atoms - [$protein_all num] - [$nucleic_all num] - [$water_all num] - [$ion_all num]}]
    puts "    Other (ligands, lipids, etc.): $other_atoms atoms ([format "%.1f" [expr {100.0 * $other_atoms / $total_atoms}]]%)"
    
    $protein_all delete
    $nucleic_all delete
    $water_all delete
    $ion_all delete
    
    puts "========================================================================"
    puts ""
}

# Quick chain info function (from previous script, enhanced)
proc list_chains {{molid "top"}} {
    if {$molid == "top"} {
        set molid [molinfo top]
    }
    
    if {$molid == -1} {
        puts "ERROR: No molecule loaded!"
        return
    }
    
    set all [atomselect $molid "all"]
    set chains [lsort -unique [$all get chain]]
    $all delete
    
    puts "\nAvailable chains: $chains"
    puts ""
    
    foreach chain $chains {
        set chain_sel [atomselect $molid "chain $chain"]
        set resids [lsort -unique [$chain_sel get resid]]
        set num_residues [llength $resids]
        
        puts "Chain $chain:"
        puts "  Residues: $num_residues"
        puts "  Range: [lindex $resids 0] to [lindex $resids end]"
        $chain_sel delete
    }
}

# Function to export structure info to file
proc export_structure_info {filename {molid "top"}} {
    if {$molid == "top"} {
        set molid [molinfo top]
    }
    
    if {$molid == -1} {
        puts "ERROR: No molecule loaded!"
        return
    }
    
    # Redirect output to file
    set output [open $filename w]
    
    # Save current stdout
    set old_stdout stdout
    
    # Temporarily redirect puts to file (this is a workaround)
    # We'll capture the analysis by running it and writing to file
    
    close $output
    
    # For now, suggest using VMD's logfile command
    puts "To export structure info to a file, you can use:"
    puts "  logfile $filename"
    puts "  analyze_structure"
    puts "  logfile off"
}

puts "\n========================================================================"
puts "VMD Structure Analysis Script Loaded"
puts "========================================================================"
puts "Available commands:"
puts "  analyze_structure \[molid\]   - Comprehensive structure analysis (default: top)"
puts "  list_chains \[molid\]         - Quick chain information (default: top)"
puts ""
puts "Examples:"
puts "  analyze_structure"
puts "  analyze_structure 0"
puts "  list_chains"
puts ""
puts "To save output to file:"
puts "  logfile structure_info.txt"
puts "  analyze_structure"
puts "  logfile off"
puts "========================================================================"
