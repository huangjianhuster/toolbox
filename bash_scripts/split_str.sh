#!/bin/bash

# Your string with a specific separator
input_string="one:two:three:four"

# Split the string using parameter expansion
IFS=":"
parts=($input_string)

# Loop through the array of parts
for part in "${parts[@]}"; do
    echo "Part: $part"
done

