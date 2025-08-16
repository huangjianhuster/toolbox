#!/bin/bash

# Example float number
float_number=12.3456789

# Using printf to format the float number
printf "Formatted Float: %6.3f\n" "$float_number"



# convert string to number
full_string="12.3456789"
bstring="${full_string:1:9}" # slicing
# Convert the substring to a number
number1=$((10#$substring))

# Another number for the calculation
number2=5

# Perform addition and store the result in a variable
result_addition=$((number1 + number2))
