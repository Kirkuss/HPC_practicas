#!/bin/bash

# Define matrix sizes
sizes=(128 256 512 1024 2048 4096)

# Output file
output_file="practica_3.txt"
echo "Matrix Multiplication Execution Times" > "$output_file"

echo "Compiling the program..."
gcc -fopenmp -O2 -o practica_3 practica_3.c

if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting."
    exit 1
fi

echo "Starting tests..."

# Loop through each matrix size
for size in "${sizes[@]}"; do
    echo "Running for matrix size ${size}x${size}..."
    echo "Matrix size: ${size}x${size}" >> "$output_file"
    for i in {1..10}; do
        ./practica_3 "$size" "$size" "$size" >> "$output_file"
    done
    echo "---------------------------------" >> "$output_file"
done

echo "Tests completed. Results saved in $output_file."

