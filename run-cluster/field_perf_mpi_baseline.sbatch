#!/bin/bash
#SBATCH -J fdtd_baseline
#SBATCH -t 02:00:00
#SBATCH -o logs/baseline_%j.out
#SBATCH -e logs/baseline_%j.err

#SBATCH -N 1
#SBATCH -n 1

module load gcc-9.5.0
module load mpi/2021.8.0

srun --cpu_bind=v,map_cpu:0 ./reference_perf_test --grid-longitudinal 1024 --grid-transverse 64
srun --cpu_bind=v,map_cpu:0 ./reference_perf_test --grid-longitudinal 128 --grid-transverse 128