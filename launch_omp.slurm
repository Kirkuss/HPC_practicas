#!/bin/bash
#SBATCH --job-name=omp_grupo10
#SBATCH -t 0-10:00
#SBATCH -p workq
#SBATCH -o logs/omp_g10.%j.out
#SBATCH -e logs/omp_g10.%j.err
#SBATCH --chdir=HPC_practicas/
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB

output_file="logs/results_omp_g10.txt"
time_file="logs/tiempos_omp_g10.txt"
sizes=(128 256 512 1024 2048 4096)


echo "Tamaño, Media, Desviacion" > "$output_file"
echo "Tamaño, Media, Desviacion" > "$time_file"

for size in "${sizes[@]}"; do
        echo "Tamaño: $size"
        tiempos=()

        for i in {1..10}; do
                echo "iteracion $i"
                tiempo=$( { /usr/bin/time -f "%e" ./omp_grupo10.exe $size $size $size > /dev/null; } 2>&1  )
                tiempos+=("$tiempo")
                echo "$size, $i, $tiempo" >> "$time_file"
        done

        avg=$(printf "%s\n" "${tiempos[@]}" | awk '{sum+=$1} END {print sum/NR}')
        std=$(printf "%s\n" "${tiempos[@]}" | awk -v mean="$avg" '{sum+=($1-mean)^2} END {print sqrt(sum/NR)}')

        echo "$size, $avg, $std" >> "$output_file"
done

#echo "Finished with job $SLURM_JOBID
