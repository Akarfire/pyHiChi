#!/bin/bash
#SBATCH -J fdtd_uniform
#SBATCH -t 02:00:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err 

module load gcc-9.5.0
module load mpi/2021.8.0

# -N and -n from arguments
srun --cpu_bind=v,map_cpu:0 ./mpi_perf_tests --grid-longitudinal 128 --grid-transverse 128 --split-mode uniform