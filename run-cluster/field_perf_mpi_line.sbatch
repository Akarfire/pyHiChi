#!/bin/bash
#SBATCH -J fdtd_line
#SBATCH -t 00:10:00
#SBATCH -o logs/perf_%j.out
#SBATCH -e logs/perf_%j.err

export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun --cpu_bind=cores ../build/unix_makefiles/src/MPI_PerformanceTests/mpi_perf_tests --grid-longitudinal 1024 --grid-transverse 64 --split-mode line