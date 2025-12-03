#!/usr/bin/env bash

# file i/o
ensure_dir() {
  local dir="$1"
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
  fi
}

# remove if already exist
remove_path() {
    local target="$1"

    if [ -e "$target" ]; then
        rm -rf "$target"
        echo "Removed: $target"
    else
        echo "Path does not exist: $target"
    fi
}

# usage: remove_suffix filepath
remove_suffix() {
    local filepath="$1"
    echo "${filepath%.*}"
}


# remove gromacs old files
clean_hash_files() {
    local dir="$1"

    if [ ! -d "$dir" ]; then
        echo "Not a directory: $dir"
        return 1
    fi

    # Find and remove files starting with '#'
    find "$dir" -maxdepth 1 -type f -name '#*' -exec rm -f {} \;

    echo "Removed files starting with '#' in: $dir"
}


# usage: build_path parentdir subdir1 subdir2 ....
build_path() {
    local base_dir="$1"
    shift

    # Remove trailing slash from base directory
    base_dir="${base_dir%/}"

    local final_path="$base_dir"

    # Process all remaining arguments
    for component in "$@"; do
        # Skip empty components
        [ -z "$component" ] && continue

        # Remove leading and trailing slashes
        component="${component#/}"
        component="${component%/}"

        # Skip if component became empty after slash removal
        [ -z "$component" ] && continue

        # Append component to path
        final_path="$final_path/$component"
    done

    echo "$final_path"
}


# parallelization
run_parallel() {
    local max_jobs="$1"
    local func_name="$2"
    local input_file="$3"

    # Validate inputs
    if [ -z "$max_jobs" ] || [ -z "$func_name" ] || [ -z "$input_file" ]; then
        echo "Usage: run_parallel <max_jobs> <function_name> <input_file>"
        return 1
    fi

    if [ ! -f "$input_file" ]; then
        echo "Error: Input file '$input_file' not found"
        return 1
    fi

    # Check if function exists
    if ! declare -f "$func_name" > /dev/null; then
        echo "Error: Function '$func_name' not found"
        return 1
    fi

    # Read all jobs into array
    local jobs=()
    while IFS= read -r line; do
        # Skip empty lines
        [ -z "$line" ] && continue
        # Skip comment lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        jobs+=("$line")
    done < "$input_file"

    local total=${#jobs[@]}
    local started=0

    if [ $total -eq 0 ]; then
        echo "Error: No jobs found in input file"
        return 1
    fi

    echo "========================================="
    echo "Parallel Execution Configuration"
    echo "========================================="
    echo "Function: $func_name"
    echo "Total jobs: $total"
    echo "Max parallel jobs: $max_jobs"
    echo "Input file: $input_file"
    echo "========================================="
    echo ""

    for job in "${jobs[@]}"; do
        # Wait if we've reached the limit
        while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
            sleep 0.5
        done

        # Parse job arguments (space-separated)
        read -r -a args <<< "$job"

        # Start the job
        $func_name "${args[@]}" &

        started=$((started + 1))
        local progress=$((started * 100 / total))
        echo "[$(date '+%H:%M:%S')] [$progress%] Started job $started/$total: ${args[*]}"
    done

    echo ""
    echo "All jobs started. Waiting for completion..."
    wait
    echo ""
    echo "========================================="
    echo "All $total jobs completed successfully!"
    echo "========================================="
}

# ### EXAMPLE of parallel calculation
# # rms functions
# rms_cal() {
#     local pdb="$1"
#     local xtc="$2"
#     local ndx="$3"
#     local out="$4"
#
#     echo "17 17" | gmx_mpi rms -s ${pdb} -f ${xtc} -n ${ndx} -fit rot+trans -dt 1 -tu ns -o ${out}_TM.xvg
#     echo "17 18" | gmx_mpi rms -s ${pdb} -f ${xtc} -n ${ndx} -fit rot+trans -dt 1 -tu ns -o ${out}_ARD.xvg
#     echo "17 19" | gmx_mpi rms -s ${pdb} -f ${xtc} -n ${ndx} -fit rot+trans -dt 1 -tu ns -o ${out}_CTER.xvg
# }
#
# ensure_dir "/home2/jianhuang/projects/thermosensation/trpv4/apo/rms/"
# outpath="/home2/jianhuang/projects/thermosensation/trpv4/apo/rms/rmsd"
# # input txt
# cat > jobs.txt << EOF
# ${eqpdb} ${xtc298Kr1} ${indexf} ${outpath}_298Kr1_
# ${eqpdb} ${xtc298Kr2} ${indexf} ${outpath}_298Kr2_
# ${eqpdb} ${xtc350Kr1} ${indexf} ${outpath}_350Kr1_
# ${eqpdb} ${xtc350Kr2} ${indexf} ${outpath}_350Kr2_
# ${eqpdb} ${xtc400Kr1} ${indexf} ${outpath}_400Kr1_
# ${eqpdb} ${xtc400Kr2} ${indexf} ${outpath}_400Kr2_
# EOF
#
# # calculation
# run_parallel 6 rms_cal jobs.txt
