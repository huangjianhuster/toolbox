#!/bin/bash

# Example: Command 1 generates multiple lines
output=$(ls)

# Loop through each file and print its name
while IFS= read -r file; do
    echo "Processing file: $file"
done <<< "$output"

